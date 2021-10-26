library(scales)
library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
library(tidyverse)
library(fmsb)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

flist= snakemake@input #list.files('/mnt/hdd/common/pol/metaGWAS/colocalization/GAraw/', 'pph_BW_', full.names=T)

funk= function(x){
d= fread(x)
d= filter(d, PP.H1.abf + PP.H2.abf + PP.H3.abf + PP.H4.abf + PP.H0.abf> 0)
fname= gsub('.txt', '', gsub('pph_', '', unlist(strsplit(x, '/'))[9]))
d= separate(d, locus, into= c('chrom', 'locus'), sep= '_')
d$sloc= d$PP.H4.abf + d$PP.H3.abf
d= select(d, PP.H4.abf, sloc, locus)

names(d)= c(fname, paste0(fname, '_sloc'), 'locus')
return(d)
}

d= lapply(flist, funk)

d= reduce(d, full_join, by = "locus")

d= arrange(d, BW_maternal_effect)

# Spider plot maternal

x= as.data.frame(matrix(d$BW_maternal_effect, ncol= nrow(d)))
x=rbind(x, as.data.frame(matrix(d$BW_maternal_effect_sloc, ncol= nrow(d))))
names(x)= d$locus

rownames(x)= c('BW maternal effect', 'BW maternal effect ')

x= rbind(rep(1,nrow(d)) , rep(0,nrow(d)) , x)


png(snakemake@output[[1]], width= 60, height= 60, res= 300, units= 'mm')
par(mar=c(0,0,0,0))

radarchart(x, axistype= 0, 
 
    #custom polygon
    pcol= c(colorBlindBlack8[4], colorBlindBlack8[2]) , pfcol= c(alpha(colorBlindBlack8[4], 0.4), alpha(colorBlindBlack8[2], 0.4)) , plwd=1, pty= 32, plty= 1,
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="#525252", caxislabels= seq(0, 1, 0.25), caxisoffset= 0.1, cglwd=0.8, calcex= 0.4,
 
    #custom labels
    vlcex= 0.43
    )

dev.off()


# Spider plot fetal

x= as.data.frame(matrix(d$BW_fetal_effect, ncol= nrow(d)))
x=rbind(x, as.data.frame(matrix(d$BW_fetal_effect_sloc, ncol= nrow(d))))
names(x)= d$locus

rownames(x)= c('BW fetal effect', 'BW fetal effect ')

x= rbind(rep(1,nrow(d)) , rep(0,nrow(d)) , x)


png(snakemake@output[[2]], width= 60, height= 60, res= 300, units= 'mm')
par(mar=c(0,0,0,0))

radarchart(x, axistype= 0,

    #custom polygon
    pcol= c(colorBlindBlack8[4], colorBlindBlack8[2]) , pfcol= c(alpha(colorBlindBlack8[4], 0.4), alpha(colorBlindBlack8[2], 0.4)) , plwd=1, pty= 32, plty= 1,
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="#525252", caxislabels= seq(0, 1, 0.25), caxisoffset= 0.1, cglwd=0.8, calcex= 0.4,

    #custom labels
    vlcex= 0.43
    )

dev.off()

