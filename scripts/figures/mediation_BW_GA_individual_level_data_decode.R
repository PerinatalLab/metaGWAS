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

as= 8
as1= 8

d= fread(snakemake@input[[1]])
d$effect= 'fetal_effect'
x= fread(snakemake@input[[2]])
x$effect= 'maternal_effect'

d= rbind(d, x)

d= filter(d, !(rsID %in% c('rs7819593', 'rs41311445')))
d$Beta2= ifelse(d$Beta1< 0, -1 * d$Beta2, d$Beta2)
d$Beta1= ifelse(d$Beta1< 0, -1 * d$Beta1, d$Beta1)

d$beta_dif= with(d, (Beta2 - Beta1) / Beta1)

mor= filter(d, effect == 'maternal_effect') %>% pull(beta_dif)
barn= filter(d, effect == 'fetal_effect') %>% pull(beta_dif)

p1= ggplot() +
geom_density( mapping=aes(x = mor, y = ..density..), fill= colorBlindBlack8[3], colour= colorBlindBlack8[3]) +
annotate('text', x= 0.35, y= 0.6, label= "Maternal", color= colorBlindBlack8[3], size= as1/ .pt, fontface = 'bold') +
annotate('text', x= 0.35, y= -1, label="Fetal", color= colorBlindBlack8[8], size= as1/ .pt, fontface = 'bold') +
geom_density(mapping= aes(x = barn, y = -..density..), fill= colorBlindBlack8[8], colour= colorBlindBlack8[8]) +
  theme_cowplot(font_size = 8) +
scale_x_continuous(expand= c(0, 0)) +
  xlab("Relative difference in effect size on \nbirth weight with or without adjusting for gestational duration") +
ylab('Density') +
geom_hline(yintercept= 0, colour= 'grey') +
theme(axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
        axis.ticks= element_line(size= 0.3))


ggsave(snakemake@output[[1]], plot= p1, width= 70, height= 70, units= 'mm', dpi= 300)



p1= ggplot(d, aes(beta_dif, group= effect, fill= effect)) +
geom_hline(yintercept= 0, colour= 'black') +
geom_density(color= NA) +
annotate('text', x=-1.5, y= 0.8, label= "Maternal", color= colorBlindBlack8[3], size= as1/ .pt, fontface = 'bold') +
annotate('text', x=1, y= 0.8, label="Fetal", color= colorBlindBlack8[8], size= as1/ .pt, fontface = 'bold') +
theme_cowplot(font_size= 8) +
#scale_colour_manual(values= alpha(colorBlindBlack8[c(8,3)], 0.5), guide= 'none') +
scale_fill_manual(values= alpha(colorBlindBlack8[c(8,3)], 0.5), guide= 'none') +
scale_x_continuous(expand= c(0, 0)) +
scale_y_continuous(expand=c(0, 0.05)) +
  xlab("Relative difference in effect size on birth weight\nwith or without adjusting for gestational duration") +
ylab('Density') +
theme(axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
        axis.ticks= element_line(size= 0.3)) +
geom_vline(xintercept= 0, linetpye= 'dashed', colour= 'grey')

ggsave(snakemake@output[[3]], plot= p1, width= 70, height= 70, units= 'mm', dpi= 300)

fwrite(d, snakemake@output[[2]], sep= '\t')
