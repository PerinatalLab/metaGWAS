library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
library(ggdendro)
library(gridExtra)
library(dendextend)
library(plyr)
library(ggtree)


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

d= fread(snakemake@input[[1]])


x= fread(snakemake@input[[2]], select= c('nearestGene', 'RSID'))

d= inner_join(d, x, by= c('rsid'= 'RSID'))

d$GENE= d$nearestGene
d$GENE= with(d, ifelse(GENE== 'CDC42', 'CDC42/ WNT4', ifelse(GENE== 'HIVEP3', 'HIVEP3/ EDN2', ifelse(GENE== 'TET3', 'TET3/ DGUOK-AS1', ifelse(GENE== 'TCEA2', 'TCEA2/ OPRL1', GENE)))))
d$nearestGene= d$GENE

d$nearestGene= with(d, ifelse(rsid== 'rs3129768', 'HLA-DQA1', ifelse(rsid== 'rs5991030', 'AGTR2', ifelse(rsid== 'rs5930554', 'RAP2C', nearestGene)))) 

d$nearestGene= with(d, ifelse(rsid== 'rs6780427', 'KCNAB1', nearestGene))
d$nearestGene= with(d, ifelse(rsid== 'rs6879092', 'EBF1', nearestGene))


d$nearestGene= gsub(' ', '', d$nearestGene)
d$nearestGene= paste0("(", d$nearestGene, ")")
d$rsid_lab= with(d, paste(rsid, nearestGene))

d$beta_PT= with(d, ifelse(beta_MT<0, -1 * beta_PT, beta_PT))
d$beta_MNT= with(d, ifelse(beta_MT<0, -1 * beta_MNT, beta_MNT))
d$beta_MT= with(d, ifelse(beta_MT<0, -1 * beta_MT, beta_MT))

d= gather(d, haplotype, beta, c('beta_MT', 'beta_MNT', 'beta_PT'))

max_beta= max(abs(d$beta))


d$haplotype= with(d, ifelse(haplotype== 'beta_MT', 'Maternal\ntransmitted', ifelse(haplotype== 'beta_MNT', 'Maternal\nnon-transmitted', 'Paternal\ntransmitted')))
d$rsid_lab= factor(d$rsid_lab, levels= unique(d$rsid_lab))


d$class_name= factor(d$class_name, levels= c("Maternal", "MF SD", "MF OD", "Fetal MatT", "Fetal"))

d= d %>% arrange(class_name, desc(probability)) %>% ungroup()
d$rsid_lab= factor(d$rsid_lab, levels= unique(d$rsid_lab))

labs <- sapply(
  strsplit(levels(d$rsid_lab), " "), 
  function(x) parse(text = paste0(x[1], "~italic('", x[2], "')"))
)

p1= ggplot(d, aes(rsid_lab, haplotype, fill= beta)) +
  theme_cowplot(8) +
  geom_tile() +
  scale_fill_gradient2(low= colorBlindBlack8[4], high= colorBlindBlack8[2], mid= 'white', limits= c(-max_beta, max_beta), guide= 'none') +
  coord_equal() +
  scale_x_discrete(labels= labs) +
  theme(axis.title= element_blank(),
        axis.ticks= element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"),
        text= element_text(size= 9/ .pt),
        axis.text.y= element_text(hjust= 0.5),
	axis.text.x= element_text(angle= 45, hjust= 1),
        axis.line = element_line(colour = 'black', size = 0.2)) +
  geom_text_repel(data= filter(d, haplotype== 'Paternal\ntransmitted'), aes(x= rsid_lab, y= 4,
                label= round(probability, 2)),  direction= 'y', size= 8/ .pt, box.padding = 0.01)  

ggsave(snakemake@output[[1]], plot= p1, width= 185, height= 60, units= 'mm', dpi= 300)

p1= ggplot(d, aes(rsid_lab, haplotype, fill= beta)) +
  theme_cowplot(8) +
  geom_tile() +
  scale_fill_gradient2(low= colorBlindBlack8[4], high= colorBlindBlack8[2], mid= 'white', limits= c(-max_beta, max_beta), name= 'Effect size') +
  coord_equal() +
scale_x_discrete(labels= labs) +  
theme(axis.title= element_blank(),
        axis.ticks= element_blank(),
        plot.margin = margin(0, 9, 0,0, "mm"),
        text= element_text(size= 9/ .pt),
        axis.text.y= element_text(hjust= 0.5),
        axis.line = element_line(colour = 'black', size = 0.2)) +
  geom_text_repel(data= filter(d, haplotype== 'Paternal\ntransmitted'), aes(x= rsid_lab, y= -0.05,
                                                                                label= round(probability, 2)), direction= "y" ,
                  size= 6.5/ .pt) 
ggsave(snakemake@output[[2]], plot= p1, width= 185, height= 100, units= 'mm', dpi= 300)

fwrite(d, snakemake@output[[3]], sep= '\t')
