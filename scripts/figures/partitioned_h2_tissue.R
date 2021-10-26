
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

d= fread('/mnt/hdd/common/pol/metaGWAS/LDscore/h2_cts/chromatin/GAraw.cell_type_results.txt')

nr= nrow(d)

d= arrange(d, Coefficient_P_value)

d= d[1:10, ]


d= separate(d, Name, into= c('Tissue', 'Mark'), sep= '__')
d$description= paste0(d$Tissue, ' (', d$Mark, ')')


d= arrange(d, desc(Coefficient_P_value))

d$description= factor(d$description, levels= unique(d$description))



p1= ggplot(data=d, aes(x= description, y= -log10(Coefficient_P_value))) +
geom_col(fill=colorBlindBlack8[2], alpha= 0.6) +
theme_cowplot(font_size= 10) +
ylab('Enrichment -log10(pvalue)') +
theme(axis.title.y=element_blank()) +
coord_flip()


d= fread('/mnt/hdd/common/pol/metaGWAS/LDscore/h2_cts/gene_expr/GAraw.cell_type_results.txt')


d$Name= gsub(".*\\.","", d$Name)

d= arrange(d, Coefficient_P_value)

d= d[1:10, ]

ggplot(data=d, aes(x= Name, y= -log10(Coefficient_P_value))) +
geom_col(fill=colorBlindBlack8[2], alpha= 0.6) +
theme_cowplot(font_size= 10) +
ylab('Enrichment -log10(pvalue)') +
theme(axis.title.y=element_blank()) +
coord_flip()


x= plot_grid(p1, p2)
