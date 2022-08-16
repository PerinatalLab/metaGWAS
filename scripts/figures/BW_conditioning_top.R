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

z= fread(snakemake@input[[1]])

z$SNP= with(z, ifelse(ref> eff, paste(chr, pos, eff, ref, sep= ':'), paste(chr, pos, ref, eff, sep= ':')))

funk= function(infile){
d= fread(infile)
names(d)[1:11]= names(d)[2:12]
d=d[, 1:11]

d$bC= ifelse(d$b< 0, -1 * d$bC, d$bC)
d$b= ifelse(d$b< 0, -1 * d$b, d$b)

d$GWAS= ifelse(grepl('BW_maternal_effect_GA', infile), 'BW_maternal_GA', 'BW_fetal_GA')

var= ifelse(grepl('BW_maternal_effect_GA', infile), 'Maternal Only', 'Fetal Only')
temp_z= z[z$origin== var, ]

d= filter(d, SNP %in% temp_z$SNP)

return(d)

}

df_list= lapply(snakemake@input[grepl('BW', snakemake@input)], funk)

d= do.call('rbind', df_list)

d$beta_dif= with(d, (bC - b) / b)


mor= filter(d, GWAS== 'BW_maternal_GA') %>% pull(beta_dif)
barn= filter(d, GWAS== 'BW_fetal_GA') %>% pull(beta_dif)

p1= ggplot() +
geom_density( mapping=aes(x = mor, y = ..density..), fill= colorBlindBlack8[3], colour= colorBlindBlack8[3]) +
annotate('text', x= 0.1, y= 3, label= "Maternal", color= colorBlindBlack8[3], size= as1/ .pt, fontface = 'bold') +
annotate('text', x= 0.1, y= -15, label="Fetal", color= colorBlindBlack8[8], size= as1/ .pt, fontface = 'bold') +
geom_density(mapping= aes(x = barn, y = -..density..), fill= colorBlindBlack8[8], colour= colorBlindBlack8[8]) +
  theme_cowplot(font_size = 8) +
scale_x_continuous(expand= c(0, 0)) +
  xlab("Relative difference in effect size on \nbirth weight after conditioning") +
ylab('Density') +
geom_hline(yintercept= 0, colour= 'grey') +
theme(axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
        axis.ticks= element_line(size= 0.3))


ggsave(snakemake@output[[1]], plot= p1, width= 70, height= 70, units= 'mm', dpi= 300)

p1= ggplot(d, aes(beta_dif, group= GWAS, fill= GWAS)) +
geom_hline(yintercept= 0, colour= 'black') +
geom_density(color= NA) +
annotate('text', x=-0.55, y= 1, label= "Maternal", color= colorBlindBlack8[3], size= as1/ .pt, fontface = 'bold') +
annotate('text', x=0.1, y= 10, label="Fetal", color= colorBlindBlack8[8], size= as1/ .pt, fontface = 'bold') +
theme_cowplot(font_size= 8) +
#scale_colour_manual(values= alpha(colorBlindBlack8[c(8,3)], 0.5), guide= 'none') +
scale_fill_manual(values= alpha(colorBlindBlack8[c(8,3)], 0.5), guide= 'none') +
scale_x_continuous(expand= c(0, 0)) +
scale_y_continuous(expand=c(0, 0.5)) +
  xlab("Relative difference in effect size on \nbirth weight after conditioning") +
ylab('Density') +
theme(axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
        axis.ticks= element_line(size= 0.3))

ggsave(snakemake@output[[3]], plot= p1, width= 70, height= 70, units= 'mm', dpi= 300)

fwrite(d, snakemake@output[[2]], sep= '\t')


