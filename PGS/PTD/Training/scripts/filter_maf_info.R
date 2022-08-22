
library(data.table)
library(dplyr)

#### Load  snplist from Plink-QC-snplist

snp <- read.table(file = snakemake@input[[1]], sep="\t", header=F)

### Load info scores for information from Reference panel

#info <- fread(file = snakemake@input[[2]], sep="\t", header=T)

#snp <- read.table('hunt-cloud/mnt/work2/chrisf/meta_gest_duration/New_PRS/Training/Data/Plink_QC_snplist/all_merged.snplist', sep="\t", header=F)

### Load info scores for information from Reference panel

#info <- read.table('hunt-cloud/mnt/archive/MOBAGENETICS/genotypes-base/aux/markerinfo/all-markerinfo.gz', sep="\t", header=T)


snp <- snp %>% rename (
  "rsid" = "V1"
)


info <- read.table(file = snakemake@input[[2]], header=T)

info <- filter(info, RefPanelAF>0.01 & RefPanelAF<0.99 & INFO>0.7)

snp= semi_join(snp, info, by= c('rsid'= 'ID'))


write.table(snp, file=snakemake@output[[1]], quote=FALSE, col.names=FALSE, row.names=FALSE)
