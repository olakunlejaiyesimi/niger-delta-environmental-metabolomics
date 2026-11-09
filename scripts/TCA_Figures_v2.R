# TCA, nucleotide and redox supplementary figures and traceable statistical tables
# Inputs: updated (2) NEG/POS matrices, updated (1) annotated statistics,
#         sample_metadata.csv, and Representative_feature_selection.csv.
# Run from the directory containing the input CSVs, or change input_dir.
# Install once: install.packages(c('ggplot2','scales'))
suppressPackageStartupMessages({library(ggplot2); library(scales)})
input_dir <- '.'
output_dir <- file.path(input_dir, 'TCA_Figures_Output')
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

find_input <- function(base, suffix=NULL) {
  names <- if (is.null(suffix)) c(paste0(base,'.csv'),paste0(base,'(1).csv')) else c(paste0(base,suffix,'.csv'),paste0(base,'.csv'))
  hits <- file.path(input_dir,names)
  hits <- hits[file.exists(hits)]
  if (!length(hits)) stop('Missing input: ',base,'; searched: ',paste(names,collapse=', '))
  hits[1]
}
read_source <- function(path, latin=FALSE) read.csv(path, check.names=FALSE, stringsAsFactors=FALSE,
                                                     fileEncoding=if(latin) 'latin1' else '')
selection <- read_source(find_input('Representative_feature_selection'))
st <- read_source(find_input('complete_v2_feature_statistics_annotated','(1)'),latin=TRUE)
st <- st[st$analysis_set=='primary',,drop=FALSE]
st$Proposed_annotation <- trimws(as.character(st$Proposed_annotation))
st$Proposed_annotation[st$Proposed_annotation %in% c('','NA','nan')] <- NA_character_
meta <- read_source(find_input('sample_metadata'))
meta <- meta[as.character(meta$primary_include) %in% c('TRUE','True','true','1') &
             meta$organism %in% c('Cassava','Plantain','Waterleaf','Flies'),,drop=FALSE]
meta$region <- ifelse(meta$origin=='Niger Delta','Niger Delta','Reference')
if(anyDuplicated(meta$sample)) stop('Duplicate sample in metadata')
raw <- list(NEG=read_source(find_input('neg_relativized_QC8_removed','(2)')),
            POS=read_source(find_input('pos_relativized_QC8_removed','(2)')))
for(mode in names(raw)) {
  if(!'row m/z' %in% names(raw[[mode]])) stop('Missing row m/z column in ',mode)
  if(anyDuplicated(raw[[mode]][['row m/z']])) stop('Duplicate feature IDs in ',mode)
  if(length(setdiff(meta$sample,names(raw[[mode]])))) stop('Missing samples in ',mode)
}
stopifnot(all(c('group','name','feature_id','mode') %in% names(selection)))
if(anyDuplicated(selection$feature_id)) stop('Duplicate representative feature IDs')
if(any(grepl('^(210\\.0972_|289\\.092)',selection$feature_id))) stop('Incorrect succinate features')
# Enforce the agreed representative features explicitly.
required <- c('citrate isomer'='191.0198_15.42_28659',
              'succinic acid'='117.0194_12.97_26188',
              'malic acid'='133.0144_11.82_24561')
for(nm in names(required)) {
  hit <- selection$feature_id[tolower(selection$name)==nm]
  if(length(hit)!=1 || hit!=required[[nm]]) stop('Wrong representative for ',nm)
}
# No silent use of the obsolete annotation column. Verify at least one corrected annotation.
for(i in seq_len(nrow(selection))) {
  f <- selection[i,]
  rows <- st[st$feature_id==f$feature_id & st$mode==f$mode,,drop=FALSE]
  if(!nrow(rows)) stop('No statistics for ',f$feature_id)
  if(!any(!is.na(rows$Proposed_annotation))) {
    # Malate at RT 11.82 was confirmed by the analyst in the updated annotation.
    # Some exported statistics versions have its Proposed_annotation column blank.
    # Do not modify statistics or infer annotations for other features.
    if (f$feature_id == '133.0144_11.82_24561' && tolower(f$name) == 'malic acid') {
      warning('Malic acid representative lacks Proposed_annotation in this source CSV; using the explicitly confirmed feature selection. Check that you supplied the latest annotated statistics file.')
    } else {
      stop('Missing corrected annotation for ', f$feature_id,
           '. Check Proposed_annotation in the updated statistics CSV.')
    }
  }
}
fmt <- function(x) if(length(x)==0 || is.na(x)) 'NA' else if(x<1e-4) formatC(x,format='e',digits=1) else formatC(x,format='f',digits=3)
organisms <- c('Cassava','Plantain','Waterleaf','Flies')
all_summary <- list(); all_samples <- list(); all_labels <- list()
for(group_name in unique(selection$group)) {
  sel <- selection[selection$group==group_name,,drop=FALSE]
  dlist <- list(); slist <- list(); llist <- list()
  for(i in seq_len(nrow(sel))) {
    f <- sel[i,]; matrix <- raw[[f$mode]]
    if(is.null(matrix)) stop('Unknown mode ',f$mode)
    idx <- match(f$feature_id,as.character(matrix[['row m/z']]))
    if(is.na(idx)) stop('Missing feature ',f$feature_id,' in ',f$mode)
    v <- suppressWarnings(as.numeric(unlist(matrix[idx,meta$sample,drop=FALSE],use.names=FALSE)))
    if(length(v)!=nrow(meta) || anyNA(v) || any(v<0)) stop('Invalid intensity: ',f$feature_id)
    d <- data.frame(group=group_name,metabolite=f$name,feature_id=f$feature_id,mode=f$mode,
                    sample=meta$sample,organism=meta$organism,region=meta$region,
                    intensity=v,signal=ifelse(v>0,'Detected','Not detected'),stringsAsFactors=FALSE)
    dlist[[i]] <- d
    for(org in organisms) {
      z <- d[d$organism==org,,drop=FALSE]; nd <- z[z$region=='Niger Delta',]; ref <- z[z$region=='Reference',]
      sr <- st[st$feature_id==f$feature_id & st$mode==f$mode & st$organism==org,,drop=FALSE]
      if(nrow(sr)>1) stop('Duplicate statistics for ',f$feature_id,' ',org)
      val <- function(col) if(nrow(sr)==1 && col %in% names(sr)) suppressWarnings(as.numeric(sr[[col]][1])) else NA_real_
      annotation <- if(nrow(sr)) sr$Proposed_annotation[1] else NA_character_
      if (is.na(annotation) && f$feature_id == '133.0144_11.82_24561' && org != 'Flies') annotation <- 'malic acid'
      row <- data.frame(group=group_name,metabolite=f$name,feature_id=f$feature_id,mode=f$mode,
                        organism=org,Proposed_annotation=annotation,
                        annotation_available=!is.na(annotation),
                        n_ND=nrow(nd),n_control=nrow(ref),n_detect_ND=sum(nd$intensity>0),
                        n_detect_control=sum(ref$intensity>0),mean_ND=mean(nd$intensity),
                        mean_control=mean(ref$intensity),
                        p_limma_two_sided=val('p_limma_two_sided'),
                        q_limma_two_sided_BH=val('q_limma_two_sided_BH'),
                        p_fisher_detection_two_sided=val('p_fisher_detection_two_sided'),
                        q_fisher_detection_BH=val('q_fisher_detection_BH'),stringsAsFactors=FALSE)
      slist[[length(slist)+1]] <- row
      label <- paste0('ND ',row$n_detect_ND,'/',row$n_ND,'  |  Ref ',row$n_detect_control,'/',row$n_control,
                      '\nlimma p ',fmt(row$p_limma_two_sided),'; q ',fmt(row$q_limma_two_sided_BH),
                      '\nFisher p ',fmt(row$p_fisher_detection_two_sided),'; q ',fmt(row$q_fisher_detection_BH))
      llist[[length(llist)+1]] <- data.frame(metabolite=f$name,organism=org,label=label)
    }
  }
  dat <- do.call(rbind,dlist); summary <- do.call(rbind,slist); labs <- do.call(rbind,llist)
  all_summary[[group_name]] <- summary; all_samples[[group_name]] <- dat; all_labels[[group_name]] <- labs
  dat$metabolite <- factor(dat$metabolite,levels=sel$name)
  labs$metabolite <- factor(labs$metabolite,levels=sel$name)
  dat$organism <- factor(dat$organism,levels=organisms)
  labs$organism <- factor(labs$organism,levels=organisms)
  dat$region <- factor(dat$region,levels=c('Reference','Niger Delta'))
  positives <- dat$intensity[dat$intensity>0]
  if(!length(positives)) stop('All intensities zero in ',group_name)
  floors <- tapply(dat$intensity,dat$metabolite,function(x) {x<-x[x>0];if(length(x))min(x)/3 else NA_real_})
  floors[is.na(floors)] <- min(positives)/3
  dat$plot_value <- ifelse(dat$intensity>0,dat$intensity,as.numeric(floors[as.character(dat$metabolite)]))
  upper <- tapply(dat$plot_value,dat$metabolite,max)
  labs$y <- as.numeric(upper[as.character(labs$metabolite)])*20
  anchor <- data.frame(metabolite=factor(sel$name,levels=sel$name),
                       organism=factor(rep(organisms[1],nrow(sel)),levels=organisms),
                       x=1.5,y=as.numeric(upper[sel$name])*22)
  set.seed(20261022)
  g <- ggplot(dat,aes(x=region,y=plot_value,color=region,shape=signal))+
    geom_point(position=position_jitter(width=.12,height=0),size=2.4,stroke=.85,alpha=.9)+
    geom_blank(data=anchor,aes(x=x,y=y),inherit.aes=FALSE)+
    geom_text(data=labs,aes(x=1.5,y=y,label=label),inherit.aes=FALSE,
              size=3.8,lineheight=1.08,vjust=1,color='black')+
    facet_grid(metabolite~organism,scales='free_y',drop=FALSE)+
    scale_y_log10(labels=scales::label_scientific(digits=1),expand=expansion(mult=c(.05,.05)))+
    scale_shape_manual(values=c('Detected'=16,'Not detected'=4))+
    scale_color_manual(values=c('Reference'='#286A8B','Niger Delta'='#B45D36'))+
    labs(title=paste0(group_name,': annotated features'),
         subtitle='Primary specimens; circles = detected, crosses = original zero intensity; p = raw, q = BH-adjusted',
         x=NULL,y='Relative intensity (log10; zero at display floor)',
         shape='Feature signal',color='Sampling group')+
    theme_bw(base_size=15)+
    theme(plot.title=element_text(size=22,face='bold'),plot.subtitle=element_text(size=13),
          strip.text=element_text(size=14,face='bold'),strip.text.y=element_text(angle=0),
          axis.text.x=element_text(angle=25,hjust=1,size=12),axis.text.y=element_text(size=11),
          legend.position='bottom',panel.grid.minor=element_blank(),
          panel.spacing=grid::unit(.7,'lines'))
  prefix <- file.path(output_dir,paste0(group_name,'_representative'))
  ggsave(paste0(prefix,'.pdf'),g,width=29,height=max(11,3.65*nrow(sel)),limitsize=FALSE)
  ggsave(paste0(prefix,'.png'),g,width=29,height=max(11,3.65*nrow(sel)),dpi=300,limitsize=FALSE)
  write.csv(summary,paste0(prefix,'_statistics.csv'),row.names=FALSE,na='')
  write.csv(dat,paste0(prefix,'_sample_intensities.csv'),row.names=FALSE)
  write.csv(labs,paste0(prefix,'_panel_labels.csv'),row.names=FALSE)
}
write.csv(do.call(rbind,all_summary),file.path(output_dir,'Representative_statistical_summary.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,all_samples),file.path(output_dir,'Representative_sample_intensities.csv'),row.names=FALSE)
write.csv(selection,file.path(output_dir,'Representative_feature_selection.csv'),row.names=FALSE)

# Comprehensive inventory: all original primary rows for selected biochemical annotations,
# including alternative chromatographic features. Never append or invent rows.
# All feature IDs are preserved and sorted together by annotation, RT, organism.
terms <- c('atp','adp','idp','udp','gdp','udp-glucose','uridine 5\u2032-diphosphate glucose',
           'citrate','citrate isomer','citric acid','isocitric acid','isocitrate','methylcitrate',
           '2-methylcitrate','succinic acid','succinate','fumarate','fumaric acid',
           'malic acid','malate','alpha-ketoglutarate','\u03b1-ketoglutarate','2-oxoglutarate',
           'oxalyl-succinic acid','oxalylsuccinate','oxaloacetate','oxaloacetic acid',
           'nad+','nadh','gssg','oxidized glutathione')
full <- st[!is.na(st$Proposed_annotation) & tolower(st$Proposed_annotation) %in% terms,,drop=FALSE]
# User removed these two from succinic acid; exclude them from succinate inventory even
# if an older file is accidentally supplied.
full <- full[!(grepl('^(210\\.0972_|289\\.092)',full$feature_id) &
                 tolower(full$Proposed_annotation)=='succinic acid'),,drop=FALSE]
full$retention_time_min <- suppressWarnings(as.numeric(vapply(strsplit(full$feature_id,'_',fixed=TRUE),function(x) if(length(x)>1)x[2] else NA_character_,character(1))))
full$selected_for_figure <- full$feature_id %in% selection$feature_id
full$figure_label <- ifelse(full$selected_for_figure,selection$name[match(full$feature_id,selection$feature_id)],'')
full <- full[order(tolower(full$Proposed_annotation),full$retention_time_min,full$feature_id,
                   match(full$organism,organisms)),,drop=FALSE]
write.csv(full,file.path(output_dir,'All_annotated_feature_statistics.csv'),row.names=FALSE,na='')
# Feature/organism rows that have no corrected annotation remain visible in a separate
# audit, but are not falsely called annotated metabolites.
selected_rows <- st[st$feature_id %in% selection$feature_id,,drop=FALSE]
selected_rows$selected_feature_name <- selection$name[match(selected_rows$feature_id,selection$feature_id)]
selected_rows <- selected_rows[order(selected_rows$selected_feature_name,
                                     match(selected_rows$organism,organisms)),,drop=FALSE]
write.csv(selected_rows,file.path(output_dir,'Selected_features_all_organisms_annotation_audit.csv'),row.names=FALSE,na='')
writeLines(c('Figures use selected feature IDs and corrected (1) statistics with (2) abundance files.',
             'The 15.42-minute citrate isomer is the plotted representative; 16.01-minute citrate remains an alternative.',
             'No merging of separate chromatographic features.',
             'Zeros are below-detection observations, displayed at a floor only for the log axis.',
             'p = unadjusted; q = BH-adjusted, separately for limma and Fisher tests.',
             'All_annotated_feature_statistics contains only nonblank Proposed_annotation rows.',
             'Selected_features_all_organisms_annotation_audit also retains rows with blank corrected annotations.',
             'No new statistical tests are calculated; source p and q are reproduced unchanged.'),
           file.path(output_dir,'README.txt'))
message('Created figures and tables in: ',normalizePath(output_dir))
