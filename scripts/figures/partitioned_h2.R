
library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
options(warn=-1)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

d= fread(snakemake@input[[1]])

nr= nrow(d)

d= filter(d, Enrichment_p< 0.05 / (nrow(d)- 1))

d$description= with(d, ifelse(Category== 'H3K27ac_HniszL2_0', 'H3K27ac', 
			ifelse(Category== 'SuperEnhancer_HniszL2_0', 'SuperEnhancer',
			ifelse(Category== 'Backgrd_Selection_StatL2_0', 'Background selection',
			ifelse(Category== 'CpG_Content_50kbL2_0', 'CpG content',
			ifelse(Category== 'BLUEPRINT_DNA_methylation_MaxCPPL2_0', 'DNA Methylation', NA))))))

d= arrange(d, desc(Enrichment_p))

d$description= factor(d$description, levels= unique(d$description))

p1= ggplot(data=d, aes(x= description, y= -log10(Enrichment_p))) +
geom_col(fill=colorBlindBlack8[2], alpha= 0.6) +
theme_cowplot(font_size= 10) +
ylab('Enrichment -log10(pvalue)') +
theme(axis.title.y=element_blank()) +
geom_hline(yintercept= -log10(0.05/ (nr -1)), linetype= 'dashed', colour= 'grey') +
coord_flip()

p2= ggplot(data=d, aes(x= description, y= Enrichment)) +
geom_col(fill=colorBlindBlack8[4], alpha= 0.6) +
theme_cowplot(font_size= 10) +
ylab('Enrichment (h2 / proportion of SNPs)') +
theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
geom_hline(yintercept= 1, linetype= 'dashed', colour= 'grey') +
coord_flip()

x= plot_grid(p1, p2)


ggsave(snakemake@output[[1]], plot= x, height= 50, width= 140, units= 'mm', dpi= 300)

fwrite(d, snakemake@output[[2]], sep= '\t')
