library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
options(warn=-1)


d= fread(snakemake@input[[1]], h= T)


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)




d= fread(snakemake@input[[1]], h= T, select= c('ID', 'pvalue', 'EAF'))
d$MAF= ifelse(d$EAF>0.5,  1 - d$EAF, d$EAF)
d= arrange(d, pvalue)
d= d[!duplicated(d$ID), ]


d= mutate(d, maf_tertiles = ntile(MAF, 3))
m1= round(max(d[d$maf_tertiles== 1, 'MAF']), 3)
m2= round(max(d[d$maf_tertiles== 2, 'MAF']), 3)


d$maf_tertiles= factor(d$maf_tertiles, levels=c("1", "2", "3"), labels=c(paste('MAF<', m1), paste(m1,'< MAF >', m2), paste('MAF>', m2)))

df= arrange(d, pvalue) %>% group_by(maf_tertiles) %>% mutate(exp1= -log10(1:length(pvalue)/length(pvalue)))

p1= ggplot(filter(df, pvalue<0.05), aes(exp1, -log10(pvalue), color= maf_tertiles)) +
  geom_point(size= 0.4) +
scale_color_manual(values= colorBlindBlack8[c(2,4,8)])+
  geom_abline(intercept = 0, slope = 1, alpha = .5) +
labs(colour="") +
theme_cowplot(font_size= 12) +
xlab('Expected (-log10(p-value))') +
ylab('Observed (-log10(p-value))') +
theme(legend.position= 'bottom') +
guides(colour = guide_legend(override.aes = list(size=3)))

ggsave(snakemake@output[[1]], plot= p1, width= 120, height= 120, units= 'mm', dpi= 300)
