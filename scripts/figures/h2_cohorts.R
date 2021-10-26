library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
options(warn=-1)


d= fread(snakemake@input[[1]], h= T)
x= fread(snakemake@input[[2]], h= T)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)
d$trait= 'Gestational\nduration'
x$trait= 'Preterm delivery'

d= rbind(d, x)

p1= ggplot(d, aes(cohort, h2, colour= cohort)) +
  geom_point() +
geom_errorbar(aes(ymin= I(h2 - 1.96*se) , ymax= (h2 + 1.96 * se)), width=.2, position=position_dodge(.9)) +
theme_cowplot(font_size= 9) +
facet_wrap(vars(trait), ncol= 1) +
scale_fill_manual(values= colorBlindBlack8[c(8,3,2,6,7, 4, 1)], guide= 'none') +
scale_colour_manual(guide= 'none', values= colorBlindBlack8[c(8,3,2,6,7, 4, 1)]) +
xlab('Cohort') +
ylab('Common SNP heritability [95% CI]') +
theme(legend.position= 'none',
	strip.background = element_blank(),
	axis.text.x= element_text(angle= 45, hjust= 1))


ggsave(snakemake@output[[1]], plot= p1, width= 60, height= 120, units= 'mm', dpi= 300)
