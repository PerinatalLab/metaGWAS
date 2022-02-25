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

x= fread(snakemake@input[[2]])

d= rbind(d, x)

names(d)= c('Name', 'no_no', 'no_yes', 'yes_no', 'yes_yes', 'candidate_gene', 'rest_genes', 'OR', 'pvalue')
d$enrichment= d$candidate_gene / d$rest_genes

d= arrange(d, desc(pvalue))

d$description= with(d, ifelse(Name== 'pli', 'Loss-of-function intolerant',
			ifelse(Name== 'dominant', 'Dominant', 'Recessive')))

d$description= factor(d$description, levels= unique(d$description))



p1= ggplot(data=d, aes(x= description, y= -log10(pvalue))) +
geom_col(fill=colorBlindBlack8[2], alpha= 0.6) +
theme_cowplot(font_size= 8) +
ylab('Enrichment -log10(pvalue)') +
theme(axis.title.y=element_blank()) +
geom_hline(yintercept= -log10(0.05/nrow(d)), linetype= 'dashed', colour= 'grey') +
coord_flip()


ggsave(snakemake@output[[1]], plot= p1, height= 35, width= 90, dpi= 300, units= 'mm')

fwrite(d, snakemake@output[[2]], sep='\t')
