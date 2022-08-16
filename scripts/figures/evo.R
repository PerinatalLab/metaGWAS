library(data.table)
library(dplyr)
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library('showtext')


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

d= fread(snakemake@input[[1]])

d$lead_snp= with(d, ifelse(lead_snp== '1:50958027', '1:50959262', ifelse(lead_snp== '9:116929327', '9:116935764', ifelse(lead_snp== '5:157896786', '5:157895049', ifelse(lead_snp== '1:22511594', '1:22462111', lead_snp)))))

x= fread(snakemake@input[[2]])

x$lead_snp= paste(x$CHR, x$POS, sep= ':')

d= inner_join(d,x, by= 'lead_snp')

d$z_score= ifelse(d$z_score> 3.5, 3.5, d$z_score)

d$nearestGene= with(d, ifelse(nearestGene== 'CDC42', 'CDC42/ WNT4', ifelse(nearestGene== 'HIVEP3', 'HIVEP3/ EDN2', ifelse(nearestGene== 'TET3', 'TET3/ DGUOK-AS1', ifelse(nearestGene== 'TCEA2', 'TCEA2/ OPRL1', nearestGene)))))

d= filter(d, !(annotation %in% c('B2', 'geva_allele_age')))

d$annotation= with(d, ifelse(annotation== 'argweave', 'ARGWEAVE', 
		ifelse(annotation== 'betascore', 'Beta score',
		ifelse(annotation== 'B2', '', 
		ifelse(annotation== 'fst_eas_afr', 'Fst AFR-EAS',
		ifelse(annotation== 'fst_eur_afr', 'Fst AFR-EUR',
		ifelse(annotation== 'fst_eur_eas', 'Fst EAS-EUR',
		ifelse(annotation== 'gerp', 'GERP',
		ifelse(annotation== 'geva_allele_age', 'Alelle age',
		ifelse(annotation== 'iES_Sabeti', 'iES',
		ifelse(annotation== 'linsigh', 'LINSIGHT',
		ifelse(annotation== 'phastCon100', 'phastCONS100',
		ifelse(annotation== 'phyloP100', 'PhyloP',
		ifelse(annotation== 'xpehh_afr2_eas', 'XPEHH AFR-EAS',
		ifelse(annotation== 'xpehh_afr2_eur', 'XPEHH AFR-EUR',
		'XPEHH EAS-EUR')))))))))))))))

p1= ggplot(d, aes(annotation, nearestGene, fill= z_score)) +
geom_tile(colour = "white", size= 1) +
theme_cowplot(font_size= 9) +
scale_fill_gradient2(low= colorBlindBlack8[2], high= colorBlindBlack8[4], mid= 'white', limits= c(-2, 4)) +
theme(axis.text.x = element_text(angle = 45, hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
scale_x_discrete(position = "top") +
geom_text(data= filter(d, pvalue.x< 0.05), aes(annotation, nearestGene, label= '*'), size= 8/ .pt) +
theme(  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks= element_blank(),
        panel.border = element_rect(colour= 'black', fill= NA, size=1),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        axis.line= element_blank(),
	axis.text.y = element_text(face = "italic")) +
coord_equal()


ggsave(snakemake@output[[1]], plot= p1, width= 140, height= 120, units= 'mm', dpi= 300)
