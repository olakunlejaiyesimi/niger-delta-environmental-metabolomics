# Lipid peroxidation section: reproducible sample-level LC-MS plots
# Place this script alongside the three CSVs named below, setwd() to that folder,
# and run source('Lipid_Peroxidation_Supplementary_Combined.R') in RStudio.
# Install once: install.packages(c('ggplot2','dplyr','readr'))
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(readr)})

resolve_file <- function(candidates) {
  found <- candidates[file.exists(candidates)]
  if (length(found)) return(found[1])
  stop('Cannot locate input file. Tried: ', paste(candidates, collapse=', '),
       '\nCurrent working directory: ', getwd(),
       '\nUse setwd("path/to/your/data/folder") before running.')
}
raw_file <- resolve_file(c('neg_relativized_QC8_removed(1).csv',
                           'neg_relativized_QC8_removed.csv'))
meta_file <- resolve_file('sample_metadata.csv')
stats_file <- resolve_file('complete_v2_feature_statistics_annotated.csv')
out_dir <- 'Lipid_Peroxidation_Supplementary_Figure_Large_Text'
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
for (f in c(raw_file,meta_file,stats_file)) if (!file.exists(f)) stop('Missing file: ',f)

# IDs are taken from the CURRENT proposed annotations. Only RT=9.43 malonic acid.
features <- data.frame(
  feature_id=c('71.0139_10.41_22230','85.0295_10.32_22042',
               '99.0452_6.51_15607','103.0037_9.43_20559',
               '141.0558_11.57_24093'),
  compound=c('Malondialdehyde (MDA)','Methylmalondialdehyde',
             'Dimethylmalondialdehyde','Malonic acid (RT 9.43 min)',
             'Putative oxoheptenoic acid'), stringsAsFactors=FALSE)

meta <- read.csv(meta_file,check.names=FALSE,stringsAsFactors=FALSE)
meta <- meta[meta$mode=='NEG' & meta$organism %in% c('Cassava','Plantain','Waterleaf','Flies') &
             as.character(meta$primary_include) %in% c('TRUE','True','true','1'),]
if(anyDuplicated(meta$sample)) stop('Duplicated NEG sample IDs in metadata')
meta$region <- ifelse(meta$origin=='Niger Delta','Niger Delta','Reference control')
raw <- read.csv(raw_file,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8')
if (!'row m/z' %in% names(raw)) stop('Missing row m/z column')
raw <- raw[as.character(raw[['row m/z']]) %in% features$feature_id,,drop=FALSE]
if(nrow(raw)!=nrow(features)) stop('Some specified features are absent from the NEG matrix')
missing_samples <- setdiff(meta$sample,names(raw))
if(length(missing_samples)) stop('Samples absent from intensity file: ',paste(missing_samples,collapse=', '))

# Long format; zeros represent nondetection, not proven absence.
long <- do.call(rbind,lapply(seq_len(nrow(raw)),function(i){
  data.frame(feature_id=as.character(raw[['row m/z']][i]),sample=meta$sample,
             intensity=suppressWarnings(as.numeric(unlist(raw[i,meta$sample],use.names=FALSE))),
             organism=meta$organism,region=meta$region,stringsAsFactors=FALSE)
}))
if(anyNA(long$intensity)) stop('Missing/non-numeric intensity values in selected features')
if(any(long$intensity<0)) stop('Negative intensities found')
long$compound <- features$compound[match(long$feature_id,features$feature_id)]
long$organism <- factor(long$organism,levels=c('Cassava','Plantain','Waterleaf','Flies'))
long$region <- factor(long$region,levels=c('Reference control','Niger Delta'))
long$compound <- factor(long$compound,levels=features$compound)
long$detected <- long$intensity>0

# Statistics shown are the EXISTING V2 results; do not re-fit or select tests for significance.
st <- read.csv(stats_file,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='latin1')
st <- st[st$analysis_set=='primary' & st$mode=='NEG' &
         st$feature_id %in% features$feature_id &
         st$organism %in% c('Cassava','Plantain','Waterleaf','Flies'),,drop=FALSE]
for (nm in c('p_limma_two_sided','q_limma_two_sided_BH',
             'p_fisher_detection_two_sided','q_fisher_detection_BH')) st[[nm]] <- as.numeric(st[[nm]])

counts <- long %>% group_by(feature_id,compound,organism,region) %>%
  summarise(n=sum(detected),total=n(),mean_intensity=mean(intensity),.groups='drop')
write.csv(counts,file.path(out_dir,'sample_detection_and_means.csv'),row.names=FALSE)
write.csv(st,file.path(out_dir,'section_existing_statistics.csv'),row.names=FALSE)
write.csv(long,file.path(out_dir,'section_sample_intensities.csv'),row.names=FALSE)

# Statistical labels: limma abundance test and Fisher detection-frequency test.
# NA means a test was not available (e.g. no detection in either group).
format_p <- function(x) {
  if (length(x)==0 || is.na(x)) return('NA')
  if (x < 0.0001) return(format(x, scientific=TRUE, digits=2))
  sprintf('%.3f', x)
}

# Construct facet-specific text from the precomputed primary statistics.
label_rows <- do.call(rbind, lapply(seq_len(nrow(features)), function(i) {
  do.call(rbind, lapply(c('Cassava','Plantain','Waterleaf','Flies'), function(org) {
    fid <- features$feature_id[i]
    dat <- long[long$feature_id == fid & as.character(long$organism) == org,,drop=FALSE]
    nd <- dat[dat$region == 'Niger Delta',,drop=FALSE]
    ct <- dat[dat$region == 'Reference control',,drop=FALSE]
    a <- st[st$feature_id == fid & st$organism == org,,drop=FALSE]
    if (nrow(a) > 1) stop('Duplicated feature/species in statistics: ',fid,' ',org)
    if (nrow(a) == 0) stop('Missing statistics: ',fid,' ',org)
    data.frame(feature_id=fid,compound=features$compound[i],organism=org,
      label=paste0('Detected: ND ',sum(nd$detected),'/',nrow(nd),
                   '; Ctrl ',sum(ct$detected),'/',nrow(ct),
                   '\nlimma p=',format_p(a$p_limma_two_sided),
                   '; q=',format_p(a$q_limma_two_sided_BH),
                   '\nFisher p=',format_p(a$p_fisher_detection_two_sided),
                   '; q=',format_p(a$q_fisher_detection_BH)),
      stringsAsFactors=FALSE)
  }))
}))
label_rows$organism <- factor(label_rows$organism,levels=levels(long$organism))
label_rows$compound <- factor(label_rows$compound,levels=levels(long$compound))
write.csv(label_rows,file.path(out_dir,'figure_statistical_labels.csv'),row.names=FALSE)

# A single display floor per feature is shared by cassava, plantain, waterleaf and flies.
# The y-axis for the four organism panels of a compound is FIXED.
# Non-detections are drawn as crosses at the floor, NOT assigned real intensity.
min_positive <- tapply(long$intensity, long$feature_id, function(v) {
  z <- v[is.finite(v) & v > 0]
  if (length(z)) min(z) / 3 else NA_real_
})
# For features with no detection anywhere, use a common harmless display floor.
fallback <- min(long$intensity[long$intensity > 0]) / 3
min_positive[is.na(min_positive)] <- fallback
long$plot_value <- ifelse(long$detected,long$intensity,
                          as.numeric(min_positive[as.character(long$feature_id)]))
long$signal <- ifelse(long$detected,'Detected','Not detected')

# One consolidated supplementary figure only. Each metabolite row shares a
# common y scale across all four species; different metabolites have separate
# scales because their ionization responses differ.
set.seed(20261022)

# Complete figure: counts, limma and Fisher statistics in every panel.
lab_overview <- label_rows
lab_overview$facet <- paste0(as.character(lab_overview$organism),'\n',
                             gsub('.*Detected: (ND [^;]+; Ctrl [^\\n]+).*','\\1',lab_overview$label))
# Five metabolite rows x four organism columns; free_y in facet_grid is
# shared across columns WITHIN each row (not independent panel scaling).
lab_overview$y <- NA_real_
for (i in seq_len(nrow(lab_overview))) {
  x <- long[long$feature_id==lab_overview$feature_id[i] &
            as.character(long$organism)==as.character(lab_overview$organism[i]),,drop=FALSE]
  lab_overview$y[i] <- max(x$plot_value)*2.5
}
p <- ggplot(long,aes(x=region,y=plot_value,color=region,shape=signal)) +
  geom_point(position=position_jitter(width=.13,height=0),size=2.5,stroke=0.85,alpha=.85) +
  geom_text(data=lab_overview,aes(x=1.5,y=y,label=label),inherit.aes=FALSE,
            size=3.6,lineheight=1.12,color='black',vjust=1) +
  facet_grid(compound~organism,scales='free_y') +
  scale_y_log10(labels=scales::label_scientific(digits=1),expand=expansion(mult=c(.08,.32))) +
  scale_shape_manual(values=c('Detected'=16,'Not detected'=4)) +
  scale_color_manual(values=c('Reference control'='#2A668A','Niger Delta'='#B15A30')) +
  labs(title='Lipid peroxidation-associated features in vegetation and flies',
       subtitle='Negative LC-MS; primary samples; fly controls include Ibadan and Port Harcourt',
       x=NULL,y='Relative intensity (log scale; zeros at display floor)',shape='Signal') +
  theme_bw(base_size=14) +
  theme(legend.position='bottom',axis.text.x=element_text(angle=30,hjust=1,size=12),
        strip.text=element_text(size=13,face='bold'),
        axis.text.y=element_text(size=11),
        plot.title=element_text(size=19,face='bold'),
        plot.subtitle=element_text(size=13),
        legend.text=element_text(size=12),legend.title=element_text(size=12),
        panel.grid.minor=element_blank())
ggsave(file.path(out_dir,'Supplementary_Figure_Lipid_Peroxidation.pdf'),p,width=23,height=22)
ggsave(file.path(out_dir,'Supplementary_Figure_Lipid_Peroxidation.png'),p,width=23,height=22,dpi=300)
writeLines(c('Statistics shown on every panel: limma two-sided raw p and BH q;',
 'Fisher exact detection raw p and BH q; ND and control detection counts.',
 'The tests are COPIED from the annotated primary V2 statistics (not refitted).',
 'This script generates ONE consolidated supplementary figure, PDF and PNG.',
 'Original zeros are crosses at a feature-specific plotting floor, not measured intensity.',
 'No pseudocount is used to change positive intensities or infer statistics.',
 'Only malonic acid at retention time 9.43 min is included.',
 'Fly controls combine two Ibadan and two Port Harcourt residential specimens.',
 'Within each compound row, all four species panels share the same log y-axis and zero display floor.',
 'The absence of a detected signal is not proof of chemical absence.',
 capture.output(sessionInfo())),file.path(out_dir,'README_and_R_session.txt'))
message('Saved consolidated supplementary figure in: ',normalizePath(out_dir))
