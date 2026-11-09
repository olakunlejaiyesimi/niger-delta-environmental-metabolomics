# Supplementary figure: glycolysis and pyruvate-associated annotated features
# Uses authoritative Proposed_annotation in complete_v2_feature_statistics_annotated.csv
# Uses primary biological specimens and the original relative-intensity matrices.
# Set working directory to folder containing both mode matrices, metadata and stats.
# Install once: install.packages(c('ggplot2','dplyr','scales'))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(scales)})

locate <- function(stem) {
  hits <- list.files('.', pattern='\\.csv$', full.names=TRUE)
  exact <- hits[basename(hits) %in% c(paste0(stem,'.csv'),paste0(stem,'(1).csv'))]
  if(length(exact)) return(exact[1])
  hits <- hits[grepl(paste0('^',stem,'(\\\\([0-9]+\\\\))?\\.csv$'),basename(hits))]
  if(length(hits)) return(hits[1])
  stop('Cannot find ',stem,'.csv in ',getwd(),'. Check getwd() and list.files().')
}
# Prefer exact names; fallback handles (1) file suffix.
find_csv <- function(base) {
  x <- c(paste0(base,'.csv'),paste0(base,'(1).csv'))
  y <- x[file.exists(x)]
  if(length(y)) return(y[1])
  stop('Missing ',base,' CSV in ',getwd(),'. Files present: ',paste(list.files(pattern='\\.csv$'),collapse=', '))
}
neg_file <- find_csv('neg_relativized_QC8_removed')
pos_file <- find_csv('pos_relativized_QC8_removed')
meta_file <- find_csv('sample_metadata')
stats_file <- find_csv('complete_v2_feature_statistics_annotated')
output_dir <- 'Glycolysis_Supplementary_Figure'
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

# Feature IDs were verified against the CURRENT proposed annotation column.
# IMPORTANT: sample_metadata.csv contains NEG-labelled rows only; metadata
# applies to both NEG and POS runs because the sample columns match.
selected <- data.frame(
 section=rep('Impact on glycolytic pathway',10),
 compound=c('Glucose','Glucose-6-phosphate','Mannose-6-phosphate (in-source)',
            'Fructose-1,6-bisphosphate','Glyceraldehyde-3-phosphate',
            '2-Phosphoglyceric acid','Pyruvic acid','Alanine',
            'Lactic acid','Putative lactic acid isomer'),
 mode=c('NEG','NEG','POS','NEG','NEG','NEG','NEG','NEG','NEG','NEG'),
 feature_id=c('179.0561_3.1_8110','259.0221_14.9_28239',
              '243.0265_14.08_31396','338.9883_16.12_28937',
              '168.9908_16.37_29027','184.9856_14.81_28150',
              '87.0088_4.26_11587','88.0404_12.55_25625',
              '89.0245_10.35_22097','89.0245_2.83_6975'),
 stringsAsFactors=FALSE)
# Verify the requested IDs against the authoritative corrected annotation.

meta <- read.csv(meta_file,check.names=FALSE,stringsAsFactors=FALSE)
meta <- meta[as.character(meta$primary_include) %in% c('TRUE','True','true','1') &
               meta$organism %in% c('Cassava','Plantain','Waterleaf','Flies'),,drop=FALSE]
meta$region <- ifelse(meta$origin=='Niger Delta','Niger Delta','Reference')
meta$organism <- as.character(meta$organism)
raws <- list(NEG=read.csv(neg_file,check.names=FALSE,stringsAsFactors=FALSE),
             POS=read.csv(pos_file,check.names=FALSE,stringsAsFactors=FALSE))
for(mode in names(raws)) {
 if(!'row m/z' %in% names(raws[[mode]])) stop('Missing row m/z column in ',mode)
 if(anyDuplicated(raws[[mode]][['row m/z']])) stop('Duplicate feature IDs in ',mode)
}
st <- read.csv(stats_file,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='latin1')
st$Proposed_annotation <- trimws(as.character(st$Proposed_annotation))
st <- st[st$analysis_set=='primary',,drop=FALSE]
for(i in seq_len(nrow(selected))) {
  match_row <- st[st$feature_id==selected$feature_id[i] & st$mode==selected$mode[i],,drop=FALSE]
  if(nrow(match_row)==0) stop('Selected feature absent from statistics: ',selected$feature_id[i])
  cat('Feature verified: ',selected$compound[i],' [',selected$mode[i],'] ',selected$feature_id[i],
      ' | Proposed_annotation: ',unique(match_row$Proposed_annotation)[1],'\n',sep='')
}


p_fmt <- function(x) {
 if(length(x)==0 || is.na(x)) return('NA')
 if(x<0.0001) return(formatC(x,format='e',digits=1))
 formatC(x,format='f',digits=3)
}

make_section <- function(section,organisms,width,height) {
 sel <- selected[selected$section==section,,drop=FALSE]
 all_long <- list(); labels <- list(); audit <- list()
 for(i in seq_len(nrow(sel))) {
   f <- sel[i,]; raw <- raws[[f$mode]]
   hit <- which(as.character(raw[['row m/z']])==f$feature_id)
   if(length(hit)!=1) stop('Feature not found exactly once: ',f$mode,' ',f$feature_id)
   # Metadata has one record per specimen (mode is NEG only), but the same
   # specimens are measured in both NEG and POS intensity matrices.
   # Do not filter metadata by mode: doing so drops ALL positive-mode samples.
   mm <- meta[meta$organism %in% organisms,,drop=FALSE]
   if(nrow(mm)==0) stop('No primary metadata specimens for ',f$compound)
   if(anyDuplicated(mm$sample)) stop('Duplicate metadata sample for ',f$mode)
   absent <- setdiff(mm$sample,names(raw))
   if(length(absent)) stop('Samples absent from ',f$mode,' matrix: ',paste(absent,collapse=', '))
   v <- suppressWarnings(as.numeric(unlist(raw[hit,mm$sample,drop=FALSE],use.names=FALSE)))
   if(length(v)!=nrow(mm) || anyNA(v) || any(v<0)) stop('Invalid intensity data for ',f$feature_id)
   df <- data.frame(section=section,compound=f$compound,mode=f$mode,feature_id=f$feature_id,
                    sample=mm$sample,organism=mm$organism,region=mm$region,intensity=v,
                    detected=v>0,stringsAsFactors=FALSE)
   all_long[[i]] <- df
   for(org in organisms) {
     sub <- df[df$organism==org,,drop=FALSE]
     if(nrow(sub)==0) stop('No samples for ',org,' / ',f$compound)
     nd <- sub[sub$region=='Niger Delta',,drop=FALSE]
     ct <- sub[sub$region=='Reference',,drop=FALSE]
     stat <- st[st$mode==f$mode & st$feature_id==f$feature_id & st$organism==org,,drop=FALSE]
     if(nrow(stat)!=1) stop('Expected one primary stats row for ',f$compound,' ',org)
     # Do not claim significance based on test selection; show both raw and adjusted.
     labels[[length(labels)+1]] <- data.frame(compound=f$compound,organism=org,
       label=paste0('ND ',sum(nd$detected),'/',nrow(nd),'  |  Ref ',sum(ct$detected),'/',nrow(ct),
       '\nlimma p ',p_fmt(as.numeric(stat$p_limma_two_sided)),'; q ',p_fmt(as.numeric(stat$q_limma_two_sided_BH)),
       '\nFisher p ',p_fmt(as.numeric(stat$p_fisher_detection_two_sided)),'; q ',p_fmt(as.numeric(stat$q_fisher_detection_BH))),
       stringsAsFactors=FALSE)
     audit[[length(audit)+1]] <- data.frame(compound=f$compound,mode=f$mode,feature_id=f$feature_id,
       organism=org,ND_detected=sum(nd$detected),ND_total=nrow(nd),
       Ref_detected=sum(ct$detected),Ref_total=nrow(ct),
       ND_mean=mean(nd$intensity),Ref_mean=mean(ct$intensity),
       limma_p=as.numeric(stat$p_limma_two_sided),limma_q=as.numeric(stat$q_limma_two_sided_BH),
       fisher_p=as.numeric(stat$p_fisher_detection_two_sided),fisher_q=as.numeric(stat$q_fisher_detection_BH))
   }
 }
 dat <- do.call(rbind,all_long); labs <- do.call(rbind,labels)
 dat$compound <- factor(dat$compound,levels=sel$compound)
 labs$compound <- factor(labs$compound,levels=sel$compound)
 dat$organism <- factor(dat$organism,levels=organisms)
 labs$organism <- factor(labs$organism,levels=organisms)
 dat$region <- factor(dat$region,levels=c('Reference','Niger Delta'))
 # A common zero-display floor and identical y limits across organisms PER compound.
 # facet_grid(scales='free_y') shares one y scale across each compound row.
 floors <- tapply(dat$intensity,dat$compound,function(x){x<-x[x>0]; if(length(x)) min(x)/3 else NA_real_})
 overall <- min(dat$intensity[dat$intensity>0])/3
 floors[is.na(floors)] <- overall
 dat$plot_value <- ifelse(dat$detected,dat$intensity,as.numeric(floors[as.character(dat$compound)]))
 dat$signal <- ifelse(dat$detected,'Detected','Not detected')
 # Place statistics in an annotation band above the observations, at the same
 # relative height for each species within a compound (avoids overlapping data).
 upper <- tapply(dat$plot_value,dat$compound,max)
 labs$y <- as.numeric(upper[as.character(labs$compound)])*20
 # Invisible anchor reserves space for all statistical annotations.
 anchor <- data.frame(compound=factor(sel$compound,levels=sel$compound),
                      organism=factor(rep(organisms[1],nrow(sel)),levels=organisms),
                      x=1.5,y=as.numeric(upper[sel$compound])*22)
 set.seed(20261022)
 g <- ggplot(dat,aes(x=region,y=plot_value,color=region,shape=signal))+
   geom_point(position=position_jitter(width=.12,height=0),size=2.4,stroke=.85,alpha=.9)+
   geom_blank(data=anchor,aes(x=x,y=y),inherit.aes=FALSE)+
   geom_text(data=labs,aes(x=1.5,y=y,label=label),inherit.aes=FALSE,
             size=3.8,lineheight=1.08,vjust=1,color='black')+
   facet_grid(compound~organism,scales='free_y',drop=FALSE)+
   scale_y_log10(labels=scales::label_scientific(digits=1),expand=expansion(mult=c(.05,.05)))+
   scale_shape_manual(values=c('Detected'=16,'Not detected'=4))+
   scale_color_manual(values=c('Reference'='#286A8B','Niger Delta'='#B45D36'))+
   labs(title='Impact on glycolytic pathway: annotated features',
        subtitle='Primary specimens; circles = detected, crosses = original zero intensity; p = raw, q = BH-adjusted',
        x=NULL,y='Relative intensity (log10; zero at display floor)',shape='Feature signal',color='Sampling group')+
   theme_bw(base_size=15)+
   theme(plot.title=element_text(size=22,face='bold'),plot.subtitle=element_text(size=13),
         strip.text=element_text(size=14,face='bold'),
         strip.text.y=element_text(angle=0),
         axis.text.x=element_text(angle=25,hjust=1,size=12),axis.text.y=element_text(size=11),
         legend.position='bottom',panel.grid.minor=element_blank(),
         panel.spacing=grid::unit(.7,'lines'))
 prefix <- 'Supplementary_Glycolysis'
 ggsave(file.path(output_dir,paste0(prefix,'.pdf')),g,width=width,height=height,limitsize=FALSE)
 ggsave(file.path(output_dir,paste0(prefix,'.png')),g,width=width,height=height,dpi=300,limitsize=FALSE)
 write.csv(dat,file.path(output_dir,paste0(prefix,'_sample_intensities.csv')),row.names=FALSE)
 write.csv(do.call(rbind,audit),file.path(output_dir,paste0(prefix,'_statistics.csv')),row.names=FALSE)
 write.csv(labs,file.path(output_dir,paste0(prefix,'_panel_labels.csv')),row.names=FALSE)
 print(g)
 invisible(g)
}

glycolysis_plot <- make_section('Impact on glycolytic pathway',c('Cassava','Plantain','Waterleaf','Flies'),width=29,height=39)
write.csv(selected,file.path(output_dir,'selected_feature_inventory.csv'),row.names=FALSE)
writeLines(c(
 'Impact on glycolytic pathway: 10 annotated features; cassava, plantain, waterleaf, flies.',
 'NEG and POS ionization matrices are processed separately. Never compare raw intensity between modes.',
 'Metadata mode is NOT used to filter samples, because both ion modes share specimen IDs.',
 'All statistics use primary sample rows from the supplied corrected statistics file.',
 'Mannose-6-phosphate (in-source) is a tentative fragment annotation, NOT confirmed free mannose-6-phosphate.',
 'Lactic acid isomer is a separate chromatographic feature, not an identified stereoisomer.',
 'Zero intensities are plotted at a floor as crosses, not imputed for statistics.',
 'Within each feature, all four organisms share the same y-axis scale.',
 'GSSG and GSH are excluded intentionally. Revisit GSSG in the oxidative stress/redox section.',
 'Do not infer glycolytic flux or oxidative stress directly from relative intensities.',
 'Vegetation references: Ibadan single site. Fly references: Ibadan (2) and Port Harcourt (2).'),
 file.path(output_dir,'README.txt'))
cat('Saved supplementary glycolysis figure to ',normalizePath(output_dir),'\n')
