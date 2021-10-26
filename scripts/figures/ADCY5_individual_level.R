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

sig_exp= filter(d, pvalue< (0.05/ 12)) %>% pull(exposure)
sig_exp_bw= filter(d, pvalue_bw< (0.05/ 12)) %>% pull(exposure)

d$genome= factor(d$genome, levels= c('Paternal', 'Fetal', 'Maternal'))

d$locus= factor(d$locus, levels= c('rs28654158 (gestational duration)', 'rs11708067 (birth weight)'), labels= c('rs28654158\n(gestational duration)', 'rs11708067\n(birth weight)'))

p1= ggplot() +
  geom_segment(data= d, aes(x=genome, xend=genome, y=beta, yend=beta_bw), color="grey") +
  geom_point(data= filter(d, !(exposure %in% sig_exp)), aes(x=genome, y=beta, size= -log10(pvalue)), shape= 21, fill= 'grey', alpha= 0.7, stroke = 1/4) +
  geom_point(data= filter(d, (exposure %in% sig_exp)), aes(x=genome, y=beta, size= -log10(pvalue)), shape= 21, fill= colorBlindBlack8[2], alpha= 0.7, stroke = 1/4) +
  geom_point(data= filter(d, !(exposure %in% sig_exp_bw)), aes(x=genome, y=beta_bw, size= -log10(pvalue_bw)), shape= 21, fill= 'grey', alpha= 0.7, stroke = 1/4) +
  geom_point(data= filter(d, (exposure %in% sig_exp_bw)), aes(x=genome, y=beta_bw, size= -log10(pvalue_bw)), shape= 21, fill= colorBlindBlack8[4], alpha= 0.7, stroke = 1/4) +
  coord_flip() + 
scale_size_continuous(range = c(.001, 3),guide= F) +
facet_wrap(vars(locus)) +
  theme_cowplot(font_size= 10) +
geom_hline(yintercept= 0, colour= 'grey', size= 0.6, linetype= 'dashed') +
theme(axis.title.y= element_blank(),
strip.background = element_blank()) +
ylab('Effect size on gestational duration, days')

ggsave(snakemake@output[[1]], p1, height= 60, width= 127, units= 'mm', dpi= 300)
