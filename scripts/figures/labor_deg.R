library(data.table)
library(dplyr)
library(cowplot)
library(ggrepel)
library('showtext')


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)



d= fread(snakemake@input[[1]])

d$Category= factor(d$Category, levels= unique(d$Category))

p1= ggplot(d, aes(Enrichment, -log10(Enrichment_p))) + 
geom_point(aes(size= Enrichment_p< 0.05/ (nrow(d)-1)), shape= 21, stroke= 0.1, fill= colorBlindBlack8[4]) +
xlab('Heritability enrichment') +
ylab('-log10(Enrichment)') +
theme_cowplot(font_size= 8) +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(panel.grid.major= element_line(colour= 'grey', size= 0.05),
legend.position = "none")  +
geom_text_repel(data= filter(d, Enrichment_p< 0.05), aes(Enrichment, -log10(Enrichment_p), label= Category), size= 8/.pt)


ggsave(snakemake@output[[1]], plot= p1, width= 120, height= 90, units= 'mm', dpi= 300)


p2= ggplot(d, aes(n_genes, -log10(Enrichment_p))) + 
geom_point(aes(size= Enrichment_p< 0.05/ (nrow(d)-1)), shape= 21, stroke= 0.1, fill= colorBlindBlack8[4]) +
xlab('Size of gene set') +
ylab('-log10(Enrichment)') +
theme_cowplot(font_size= 8) +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(panel.grid.major= element_line(colour= 'grey', size= 0.05),
legend.position = "none") +
geom_text_repel(data= filter(d, Enrichment_p< 0.05), aes(n_genes, -log10(Enrichment_p), label= Category), size= 8/.pt)
 
ggsave(snakemake@output[[2]], plot= p2, width= 90, height= 90, units= 'mm', dpi= 300)

p3= ggplot(d, aes(n_genes, -log10(Enrichment_p))) + 
geom_point(aes(size= Enrichment_p< 0.05/ (nrow(d)-1)), shape= 21, stroke= 0.1, fill= colorBlindBlack8[4]) +
xlab('Size of gene set') +
ylab('-log10(Enrichment)') +
theme_cowplot(font_size= 8) +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(panel.grid.major= element_line(colour= 'grey', size= 0.05)) +
geom_text_repel(data= filter(d, Enrichment_p< 0.05), aes(n_genes, -log10(Enrichment_p), label= Category), size= 8/.pt)
 
ggsave(snakemake@output[[3]], plot= p3, width= 90, height= 90, units= 'mm', dpi= 300)
