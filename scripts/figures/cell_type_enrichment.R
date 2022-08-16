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
x= fread(snakemake@input[[2]])
d= inner_join(d, x, by= 'Name')


d= d[sample(nrow(d)),]

d= d[order(d$Category, decreasing= F), ]

d$Name= factor(d$Name, levels= unique(d$Name))

d$Name2= gsub('_', ' ', gsub("^.*\\.","", d$Name))
d$Name2= factor(d$Name2, levels= unique(d$Name2))

p1= ggplot(d, aes(Name2, -log10(Coefficient_P_value), colour= Category, fill= Category)) + 
geom_point(size= 2, shape= 21, stroke= 0.1) +
xlab('Tissues') +
ylab('-log10(Enrichment)') +
theme_cowplot(font_size= 8) +
geom_hline(yintercept= -log10(0.05), colour= '#d9d9d9') +
theme(axis.text.x = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.05),
	panel.grid.major.x= element_blank(),
	legend.position="none") +
geom_text_repel(data= filter(d, Coefficient_P_value< 0.05), aes(Name2, -log10(Coefficient_P_value), colour= Category, label= Name2, show_guide = FALSE))


ggsave(snakemake@output[[1]], plot= p1, width= 120, height= 90, units= 'mm', dpi= 300)

p2= ggplot(d, aes(Name2, -log10(Coefficient_P_value), colour= Category, fill= Category)) +
geom_point(size= 2, shape= 21, stroke= 0.1) +
xlab('Tissues') +
ylab('-log10(Enrichment)') +
theme_cowplot(font_size= 8) +
geom_hline(yintercept= -log10(0.05), colour= '#d9d9d9') +
theme(axis.text.x = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.05),
        panel.grid.major.x= element_blank()) +
geom_text_repel(data= filter(d, Coefficient_P_value< 0.05), aes(Name2, -log10(Coefficient_P_value), colour= Category, label= Name2), show_guide = FALSE)

ggsave(snakemake@output[[2]], plot= p2, width= 120, height= 90, units= 'mm', dpi= 300)
