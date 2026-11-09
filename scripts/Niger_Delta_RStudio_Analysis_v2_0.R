# Niger Delta LC-MS differential abundance analysis — Version 2.0
# RStudio: set working directory to folder containing this script and four input CSVs.
# Required: limma (Bioconductor), ggplot2 (CRAN).
# Install once:
# install.packages("ggplot2")
# if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install("limma")
# Run: source("Niger_Delta_RStudio_Analysis_v2_0.R")
# NOTE: All tests are exploratory at the individual-sample level; not site-adjusted.
# Two-sided moderated limma is PRIMARY; directional one-sided tests SUPPLEMENTARY.
# Biological absence cannot be inferred from zero relative intensity.
suppressPackageStartupMessages({
  if(!requireNamespace("limma",quietly=TRUE)) stop("Install limma with BiocManager::install('limma')")
  if(!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2 with install.packages('ggplot2')")
})
out <- "Niger_Delta_v2_results"
dir.create(out,showWarnings=FALSE)
dir.create(file.path(out,"figures"),showWarnings=FALSE)
inputs <- c(NEG="neg_relativized_QC8_removed.csv",
            POS="pos_relativized_QC8_removed.csv")
meta_file <- "sample_metadata.csv"
ann_file <- "highlighted_small_molecules_final_OJ_10.12.2026.csv"
stopifnot(all(file.exists(c(inputs,meta_file))))
meta <- read.csv(meta_file,stringsAsFactors=FALSE,check.names=FALSE)
required <- c("sample","organism","origin","primary_include")
if(!all(required %in% names(meta))) stop("Metadata missing columns: ",paste(setdiff(required,names(meta)),collapse=", "))
meta$primary_include <- tolower(as.character(meta$primary_include)) %in% c("true","t","1","yes")
meta$organism <- trimws(meta$organism)
meta$origin <- trimws(meta$origin)
if(anyDuplicated(meta$sample)) stop("Duplicate sample IDs in metadata.")
if(!all(meta$organism %in% c("Cassava","Plantain","Waterleaf","Flies"))) stop("Unexpected organism names in metadata")
if(!all(meta$origin %in% c("Niger Delta","Ibadan","Port Harcourt residential"))) stop("Unexpected origin labels: ",paste(unique(meta$origin),collapse=", "))
if(any(!meta$primary_include & !meta$sample %in% c("vdc111-113","vdn01WL","vdn01PL")))
  stop("Unexpected primary exclusion; check metadata.")
if(!all(c("vdc111-113","vdn01WL","vdn01PL") %in% meta$sample[!meta$primary_include]))
  stop("Expected pooled and two unknown-site samples excluded in primary metadata.")
if(file.exists(ann_file)) {
  aa <- read.csv(ann_file,stringsAsFactors=FALSE,check.names=FALSE)
  ann <- data.frame(feature_id=as.character(aa[["row m/z"]]),
                    annotation=as.character(aa[[1]]),stringsAsFactors=FALSE)
  ann <- ann[!is.na(ann$feature_id) & ann$feature_id!="Label" & nzchar(ann$feature_id) &
             !is.na(ann$annotation) & nzchar(trimws(ann$annotation)),]
  ann <- ann[!duplicated(ann$feature_id),]
} else {
  ann <- data.frame(feature_id=character(),annotation=character())
  warning("Annotation file absent; annotation column will be blank.")
}
# Transform: log2(1+x) of supplied relative intensities. This avoids enormous
# log fold changes for exclusive features, but coefficients are differences on
# log2(1+x) scale, NOT log2 fold changes. See methods note.
# Filtering: feature detected in >= 2 samples in EITHER group (predefined,
# symmetric), and finite intensity for all samples. Detection is x>0.
# For group-specific detections, report Fisher test separately.
# No site-level or causal inference from p-values.
result_list <- list(); count_list <- list(); robustness <- list()
for(mode in names(inputs)) {
  dat <- read.csv(inputs[[mode]],check.names=FALSE,stringsAsFactors=FALSE)
  if(!"row m/z" %in% names(dat) || as.character(dat[["row m/z"]][1])!="Label")
    stop("Unexpected feature matrix format for ",mode)
  available <- intersect(meta$sample,names(dat))
  if(length(available)!=nrow(meta)) stop("Missing metadata samples in ",mode,": ",paste(setdiff(meta$sample,available),collapse=", "))
  mat <- as.matrix(dat[-1,available,drop=FALSE])
  suppressWarnings(storage.mode(mat)<-"numeric")
  if(any(!is.finite(mat),na.rm=TRUE)) stop("Nonfinite values found in ",mode)
  if(any(mat<0,na.rm=TRUE)) stop("Negative intensities found in ",mode)
  ids <- as.character(dat[["row m/z"]][-1])
  if(anyDuplicated(ids)) stop("Duplicate feature IDs in ",mode)
  rownames(mat) <- ids
  for(org in c("Cassava","Plantain","Waterleaf","Flies")) {
    for(which_set in c("primary","sensitivity")) {
      m <- meta[meta$organism==org & (which_set=="sensitivity" | meta$primary_include),,drop=FALSE]
      group <- ifelse(m$origin=="Niger Delta","ND","Control")
      if(length(unique(group))!=2) stop("Both groups required: ",mode," ",org)
      x <- mat[,m$sample,drop=FALSE]
      detected <- x>0
      nd <- group=="ND"; ct <- !nd
      n_nd <- sum(nd);n_ct<-sum(ct)
      n_det_nd <- rowSums(detected[,nd,drop=FALSE],na.rm=TRUE)
      n_det_ct <- rowSums(detected[,ct,drop=FALSE],na.rm=TRUE)
      mean_nd <- rowMeans(x[,nd,drop=FALSE],na.rm=TRUE)
      mean_ct <- rowMeans(x[,ct,drop=FALSE],na.rm=TRUE)
      median_nd <- apply(x[,nd,drop=FALSE],1,median,na.rm=TRUE)
      median_ct <- apply(x[,ct,drop=FALSE],1,median,na.rm=TRUE)
      classification <- ifelse(n_det_nd==0 & n_det_ct==0,"Neither group",
                         ifelse(n_det_nd>0 & n_det_ct==0,"ND only detected",
                         ifelse(n_det_nd==0 & n_det_ct>0,"Control only detected",
                         ifelse(mean_nd>mean_ct,"ND higher",
                         ifelse(mean_nd<mean_ct,"Control higher","Equal mean")))))
      # Conservative symmetric detection filter, not based on direction.
      eligible <- (n_det_nd>=2 | n_det_ct>=2) & rowSums(is.na(x))==0
      p_two <- p_up <- p_down <- rep(NA_real_,nrow(x))
      coefficient <- rep(NA_real_,nrow(x))
      if(sum(eligible)>0) {
        design <- model.matrix(~0+factor(group,levels=c("Control","ND")))
        colnames(design)<-c("Control","ND")
        fit <- limma::lmFit(log2(1+x[eligible,,drop=FALSE]),design)
        contrast <- limma::makeContrasts(ND-Control,levels=design)
        fit <- limma::contrasts.fit(fit,contrast)
        fit <- limma::eBayes(fit,robust=FALSE)
        tt <- fit$t[,1]
        df <- fit$df.total
        p_two[eligible] <- fit$p.value[,1]
        p_up[eligible] <- stats::pt(tt,df=df,lower.tail=FALSE)
        p_down[eligible] <- stats::pt(tt,df=df,lower.tail=TRUE)
        coefficient[eligible] <- fit$coefficients[,1]
      }
      q_two <- stats::p.adjust(p_two,method="BH")
      # Directional tests are adjusted separately, and are SUPPLEMENTARY.
      q_up <- stats::p.adjust(p_up,method="BH")
      q_down <- stats::p.adjust(p_down,method="BH")
      # Fisher detection test: no p when neither group detects a feature.
      p_fisher <- rep(NA_real_,nrow(x))
      for(i in which((n_det_nd+n_det_ct)>0)) {
        tab <- matrix(c(n_det_nd[i],n_nd-n_det_nd[i],
                        n_det_ct[i],n_ct-n_det_ct[i]),nrow=2,byrow=TRUE)
        p_fisher[i] <- stats::fisher.test(tab,alternative="two.sided")$p.value
      }
      q_fisher <- stats::p.adjust(p_fisher,"BH")
      r <- data.frame(feature_id=ids,mode=mode,organism=org,analysis_set=which_set,
                      n_ND=n_nd,n_control=n_ct,n_detect_ND=n_det_nd,n_detect_control=n_det_ct,
                      fraction_detect_ND=n_det_nd/n_nd,fraction_detect_control=n_det_ct/n_ct,
                      mean_ND=mean_nd,mean_control=mean_ct,median_ND=median_nd,
                      median_control=median_ct,direction=classification,
                      eligible_limma=eligible,
                      limma_log2_1plus_difference=coefficient,
                      p_limma_two_sided=p_two,q_limma_two_sided_BH=q_two,
                      p_limma_ND_higher_one_sided=p_up,q_limma_ND_higher_BH=q_up,
                      p_limma_ND_lower_one_sided=p_down,q_limma_ND_lower_BH=q_down,
                      p_fisher_detection_two_sided=p_fisher,
                      q_fisher_detection_BH=q_fisher,
                      stringsAsFactors=FALSE)
      r$annotation <- ann$annotation[match(r$feature_id,ann$feature_id)]
      result_list[[paste(mode,org,which_set)]]<-r
      count_list[[paste(mode,org,which_set)]]<-data.frame(
        mode=mode,organism=org,analysis_set=which_set,
        n_ND=n_nd,n_control=n_ct,n_features=nrow(r),n_limma_tested=sum(eligible),
        ND_higher_descriptive=sum(classification=="ND higher"),
        control_higher_descriptive=sum(classification=="Control higher"),
        ND_only_detected=sum(classification=="ND only detected"),
        control_only_detected=sum(classification=="Control only detected"),
        limma_BH_significant_ND_elevated=sum(q_two<.05 & coefficient>0,na.rm=TRUE),
        limma_BH_significant_ND_reduced=sum(q_two<.05 & coefficient<0,na.rm=TRUE),
        one_sided_BH_significant_ND_elevated=sum(q_up<.05 & coefficient>0,na.rm=TRUE),
        one_sided_BH_significant_ND_reduced=sum(q_down<.05 & coefficient<0,na.rm=TRUE),
        fisher_BH_significant=sum(q_fisher<.05,na.rm=TRUE))
      # Volcano figure: x is transformed mean difference, not fold change.
      plotdata <- r[is.finite(r$limma_log2_1plus_difference) &
                    is.finite(r$p_limma_two_sided),]
      if(nrow(plotdata)>0) {
        plotdata$status <- ifelse(plotdata$q_limma_two_sided_BH<.05 &
                                   plotdata$limma_log2_1plus_difference>0,"ND elevated",
                            ifelse(plotdata$q_limma_two_sided_BH<.05 &
                                     plotdata$limma_log2_1plus_difference<0,"ND reduced","Not FDR significant"))
        plotdata$neglog10p <- -log10(pmax(plotdata$p_limma_two_sided,.Machine$double.xmin))
        g <- ggplot2::ggplot(plotdata,ggplot2::aes(x=limma_log2_1plus_difference,y=neglog10p,color=status))+
          ggplot2::geom_point(alpha=.55,size=.8)+
          ggplot2::scale_color_manual(values=c("ND elevated"="#B64B42","ND reduced"="#377BA8","Not FDR significant"="grey70"))+
          ggplot2::theme_bw(base_size=11)+
          ggplot2::labs(title=paste(org,mode,which_set),
                        subtitle="Exploratory sample-level limma; not site-adjusted",
                        x="Difference in mean log2(1 + relative intensity): ND - control",
                        y="-log10(two-sided raw p-value)",color=NULL)
        ggplot2::ggsave(file.path(out,"figures",paste0("volcano_",mode,"_",org,"_",which_set,".pdf")),
                        g,width=8,height=5.5)
      }
    }
  }
}
all <- do.call(rbind,result_list)
counts <- do.call(rbind,count_list)
write.csv(all,file.path(out,"complete_v2_feature_statistics.csv"),row.names=FALSE,na="")
write.csv(counts,file.path(out,"v2_summary_counts.csv"),row.names=FALSE)
write.csv(all[all$analysis_set=="primary" &
                all$limma_log2_1plus_difference>0 &
                !is.na(all$limma_log2_1plus_difference),],
          file.path(out,"ND_elevated_primary.csv"),row.names=FALSE,na="")
write.csv(all[all$analysis_set=="primary" &
                all$limma_log2_1plus_difference<0 &
                !is.na(all$limma_log2_1plus_difference),],
          file.path(out,"ND_reduced_primary.csv"),row.names=FALSE,na="")
write.csv(all[all$analysis_set=="primary" & all$direction %in%
                c("ND only detected","Control only detected"),],
          file.path(out,"group_specific_detection_primary.csv"),row.names=FALSE,na="")
primary <- all[all$analysis_set=="primary",]
sensitivity <- all[all$analysis_set=="sensitivity",]
keys <- c("feature_id","mode","organism")
rob <- merge(primary[,c(keys,"limma_log2_1plus_difference","q_limma_two_sided_BH")],
             sensitivity[,c(keys,"limma_log2_1plus_difference","q_limma_two_sided_BH")],
             by=keys,suffixes=c("_primary","_sensitivity"),all=TRUE)
write.csv(rob,file.path(out,"primary_sensitivity_comparison.csv"),row.names=FALSE,na="")
writeLines(c(
 "Niger Delta LC-MS Version 2.0 — statistical interpretation",
 "Primary: two-sided limma empirical Bayes on log2(1 + relative intensity).",
 "Directional one-sided limma and Fisher detection tests are supplementary/exploratory.",
 "limma_log2_1plus_difference is NOT a log2 fold change.",
 "BH correction is applied within each organism x ionization mode x analysis set.",
 "Direction is based on original group means; inferential direction is based on limma coefficient.",
 "Feature eligibility: >=2 detected samples in either group; all intensities finite.",
 "Detection: original recorded intensity >0. Zeros are not proof of biological absence.",
 "Pooled and unknown-site vegetation samples excluded from primary; included in sensitivity.",
 "Do not interpret these p-values as site-adjusted or causal evidence.",
 "Before publication, inspect QC, normalization, and abundance distribution assumptions.",
 "Annotation is matched by feature ID and must be checked against mode and MS/MS evidence."
 ),file.path(out,"README_methods_and_limitations.txt"))
capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
message("Completed. Outputs in: ",normalizePath(out))
