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

d= filter(d, MarkerName!= '6:32595083:G:T')

top= fread(snakemake@input[[2]])
ids= pull(top, ID)
ids= c('3:156697097:A:G', '5:158058432:G:T', ids)
d$ID= d$MarkerName

d= filter(d, ID %in% ids)


d= separate(d, MarkerName, into= c('CHR', 'POS', 'REF', 'EFF'), sep= ':')
d$beta_h1= with(d, ifelse(REF > EFF, -1 * beta_h1, beta_h1))
d$beta_h2= with(d, ifelse(REF > EFF, -1 * beta_h2, beta_h2))
d$beta_h3= with(d, ifelse(REF > EFF, -1 * beta_h3, beta_h3))
d$beta_h1= with(d, ifelse(d$beta_h2<0, d$beta_h1 * -1, d$beta_h1))
d$beta_h3= with(d, ifelse(d$beta_h2<0, d$beta_h3 * -1, d$beta_h3))
d$beta_h2= with(d, ifelse(d$beta_h2<0, d$beta_h2 * -1, d$beta_h2))

d$ID= with(d, ifelse(REF> EFF, paste(CHR, POS, EFF, REF, sep= ':'), paste(CHR, POS, REF, EFF, sep= ':')))

d= left_join(d, top[, c('ID', 'nearestGene', 'BETA', 'RSID', 'pvalue')], by= 'ID')
d$nearestGene= ifelse(is.na(d$nearestGene), 'EBF1', d$nearestGene)

d$GENE= d$nearestGene
d$GENE= with(d, ifelse(GENE== 'CDC42', 'CDC42/ WNT4', ifelse(GENE== 'HIVEP3', 'HIVEP3/ EDN2', ifelse(GENE== 'TET3', 'TET3/ DGUOK-AS1', ifelse(GENE== 'TCEA2', 'TCEA2/ OPRL1', GENE)))))
d$nearestGene= d$GENE

d$RSID= ifelse(d$ID== '5:158058432:G:T', 'rs6879092', d$RSID)
d$pvalue= ifelse(d$ID== '5:158058432:G:T', 1e-7, d$pvalue)

d$RSID= ifelse(d$ID== '3:156697097:A:G', 'rs6780427', d$RSID)
d$pvalue= ifelse(d$ID== '3:156697097:A:G', 2.633e-08, d$pvalue)
d$nearestGene= ifelse(d$ID== '3:156697097:A:G', 'KCNAB1', d$nearestGene)


d$rsid_lab= with(d, paste0(RSID, ' (', nearestGene, ')'))


ds= dist((d[, c('beta_h2', 'beta_h1', 'beta_h3')]))

hc= hclust(ds)
hc$labels= d$rsid_lab
clust= cutree(hc, k= 4)

d$clust= clust

circ1 <- ggtree(hc, layout = "circular")
df1= as.data.frame(select(d, beta_h2, beta_h1, beta_h3))
#names(df1)= c('MnT', 'MT', 'PT')
names(df1)= c('Maternal\nnon-transmitted', 'Maternal\ntransmitted', 'Paternal\ntransmitted')

rownames(df1) <- d$rsid_lab
print('
circ1= ggtree(hc) +
theme(legend.position="none") +
geom_hilight(node=21, fill= colorBlindBlack8[7], alpha=0.5) +
geom_cladelabel(21, "Mt", offset= 0.03, barsize=0, angle=90, offset.text= -1/100, hjust=0.5, fontsize=8/.pt) +
geom_hilight(node=24, fill= colorBlindBlack8[3], alpha=0.5) +
geom_cladelabel(24, "Opposite", offset=0.03, barsize=0, angle=90, offset.text=-1/100, hjust=0.5, fontsize=8/.pt) +
geom_hilight(node=27, fill=colorBlindBlack8[8], alpha=0.5) +
geom_cladelabel(27, "Maternal", offset=0.03, barsize=0, angle=90, offset.text=-1/100, hjust=0.5, fontsize=8/.pt) +
geom_hilight(node=26, fill="grey", alpha=0.5)  +
geom_cladelabel(26, "Unclassified", offset=0.03, barsize=0, angle=90, offset.text=-1/100, hjust=0.5, fontsize=8/.pt) +
theme(plot.margin=margin(0, 0, 0, 0, unit= "cm"))+
scale_x_continuous(expand= c(0, NA), limits= c(0, 3.1))

grp= split(d$rsid_lab, d$clust)
max_beta= max(abs(as.numeric(stack(d)$values)))
 

p1= groupOTU(circ1, grp, "RSID") + 
aes(color= RSID ) + 
theme(legend.position="none", 
plot.margin=margin(0, 0, 0, 0, unit= "cm")) +
scale_colour_manual(values=c(colorBlindBlack8[c(1, 8, 3)], "grey", colorBlindBlack8[7]), guide= "none")


ctree= gheatmap(p1, df1, offset= 0, width=0.7, 
        colnames=TRUE, legend_title="genotype", colnames_angle= 0, font.size= 8/.pt, hjust= 0.5) +
    geom_tiplab(align = TRUE, linesize=0, offset= 0.9, size= 8/.pt) +
    scale_fill_gradient2(low= colorBlindBlack8[2], high= colorBlindBlack8[4], mid= "white", guide= F, limits= c(-max_beta, max_beta)) +
coord_cartesian(clip = "off") + 
  theme(plot.margin=margin(0, 0, 0, 0, unit= "cm"))

ggsave(snakemake@output[[1]], plot= ctree, width= 100, height= 120, units= "mm", dpi= 300)
')


circ1= ggtree(hc) +
theme(legend.position='none') +
geom_hilight(node=21, fill= colorBlindBlack8[7], alpha=0.5) +
#geom_cladelabel(21, "Mt", offset= 0.0, barsize=0, angle=0, offset.text= -0.15, hjust=0.5, fontsize=8/.pt) +
geom_hilight(node=24, fill= colorBlindBlack8[3], alpha=0.5) +
#geom_cladelabel(24, "Opp.", offset=0.0, barsize=0, angle=0, offset.text=-0.3, hjust=0.5, fontsize=8/.pt) +
geom_hilight(node=14, fill=colorBlindBlack8[8], alpha=0.5) +
#geom_cladelabel(27, "Maternal", offset=0.0, barsize=0, angle=0, offset.text=-0.45, hjust=0.5, fontsize=8/.pt) +
geom_hilight(node=26, fill="grey", alpha=0.5)  +
#geom_cladelabel(26, "Unclassified", offset=0.0, barsize=0, angle=0, offset.text=-0.3, hjust=0.5, fontsize=8/.pt) +
theme(plot.margin=margin(0, 0, 0, 0, unit= 'cm'))+
scale_x_continuous(expand= c(0, NA), limits= c(0, 3.1))

grp= split(d$rsid_lab, d$clust)
max_beta= max(abs(as.numeric(stack(d)$values)))


p1= groupOTU(circ1, grp, 'RSID') +
aes(color= RSID, angle= 45, hjust= 1) +
theme(legend.position="none",
plot.margin=margin(0, 0, 0, 0, unit= 'cm')) +
scale_colour_manual(values=c(colorBlindBlack8[c(1, 8, 3)], 'grey', colorBlindBlack8[7]), guide= 'none')


ctree= gheatmap(p1, df1, offset= 0, width= 0.7,
        colnames= TRUE, legend_title= "genotype", colnames_offset_y= -1, colnames_angle= 0, font.size= 8/.pt, hjust= 0.5) +
    geom_tiplab(align = FALSE, linesize= 0, offset= -0.9, size= 8/.pt, hjust= 1) +
    scale_fill_gradient2(low= colorBlindBlack8[2], high= colorBlindBlack8[4], mid= 'white', guide= F, limits= c(-max_beta, max_beta)) +
coord_cartesian(clip = 'off') +
  theme(plot.margin=margin(0, 0, 0.6, 0.5, unit= 'cm')) +
  scale_x_reverse(limits= c(2.5, 0))+ 
coord_flip(clip= 'off')
  
ggsave(snakemake@output[[1]], plot= ctree, width= 185, height= 100, units= 'mm', dpi= 300)


fwrite(d, snakemake@output[[2]], sep= '\t')
