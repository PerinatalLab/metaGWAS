library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(tidyr)
library(showtext)
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')
as= 8
as1= 9
showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)


d= fread(snakemake@input[[1]])
names(d)= c('CHR', 'POS', 'FST_EUR_AFR')

d1= fread(snakemake@input[[2]])
names(d1)= c('CHR', 'POS', 'FST_EUR_EAS')

d2= fread(snakemake@input[[3]])
names(d2)= c('CHR', 'POS', 'FST_AFR_EAS')

d= inner_join(d, d1, by= c('CHR', 'POS')) %>% inner_join(., d2, by= c('CHR', 'POS'))

d$v_ids= paste(d$CHR, d$POS, sep= ':')

z= fread(snakemake@input[[4]])

zl= gather(z, control_set, v_ids, Set_1:Set_10000)

bw_pos= c(123065778)
ga_pos= c(123112292)

zl= inner_join(zl, d[, c('v_ids', 'FST_EUR_AFR', 'FST_EUR_EAS', 'FST_AFR_EAS')], by= 'v_ids')

zl= filter(zl, Input_SNP== '3:123065778' | Input_SNP== '3:123112292')

zl$haplotype= with(zl, ifelse(Input_SNP== '3:123065778', 'Birth weight', 'Gestational duration'))

zl= zl[!duplicated(zl$v_ids), ]

zl= data.frame(zl)
d= data.frame(d)
df_list= list()
r_num= 1
for (i in c('FST_EUR_AFR', 'FST_AFR_EAS', 'FST_EUR_EAS')){

ga_pvalue=wilcox.test(zl[zl$haplotype== 'Gestational duration', i], mu= d[d$v_ids== '3:123112292', i], alternative= 'less')$p.value

m1= d[d$v_ids== '3:123112292', i]
mc1= mean(zl[zl$haplotype== 'Gestational duration', i], na.rm=T)
medc1= median(zl[zl$haplotype== 'Gestational duration', i], na.rm=T)
prop_above= prop.table(table(d[d$v_ids== '3:123112292', i]> zl[zl$haplotype== 'Gestational duration', i]))[2]
temp_df= data.frame(haplotype= 'Gestational duration', ancestries= i, FST= m1, FST_mean_controls= mc1, FST_median_controls= medc1, pvalue= ga_pvalue)


ga_pvalue= wilcox.test(zl[zl$haplotype== 'Birth weight', i], mu= d[d$v_ids== '3:123065778', i], alternative= 'less')$p.value

medc1= median(zl[zl$haplotype== 'Birth weight', i], na.rm=T)
m1= d[d$v_ids== '3:123065778', i]
mc1= mean(zl[zl$haplotype== 'Birth weight', i], na.rm=T)

temp_df2= data.frame(haplotype= 'Birth weight', ancestries= i, FST= m1, FST_mean_controls= mc1, FST_median_controls= medc1, pvalue= ga_pvalue)
temp_df= rbind(temp_df, temp_df2)
df_list[[r_num]]= temp_df

r_num= r_num + 1
}

xp= do.call('rbind', df_list)

xp$enrichment= with(xp, FST / FST_median_controls)

bw= filter(zl, haplotype== 'Birth weight') %>% select(FST_EUR_AFR, FST_EUR_EAS, FST_AFR_EAS)
ga= filter(zl, haplotype== 'Gestational duration') %>% select(FST_EUR_AFR, FST_EUR_EAS, FST_AFR_EAS)

names(bw)= c('FST_EUR_AFR_bw', 'FST_EUR_EAS_bw', 'FST_AFR_EAS_bw')

df1= cbind(bw, ga)

ga_fst= d[d$v_ids== '3:123112292', 'FST_EUR_AFR']
bw_fst= d[d$v_ids== '3:123065778', 'FST_EUR_AFR']
ga_fst_pvalue= xp[xp$haplotype== 'Gestational duration' & xp$ancestries== 'FST_EUR_AFR', 'enrichment']
bw_fst_pvalue= xp[xp$haplotype== 'Birth weight' & xp$ancestries== 'FST_EUR_AFR', 'enrichment']

p1= ggplot(df1, aes(x=x) ) +
  geom_density( aes(x = FST_EUR_AFR, y = ..density..), fill= colorBlindBlack8[4], colour= colorBlindBlack8[4]) +
annotate('text', x=0.6, y= 10, label="Gestational \nduration", color= colorBlindBlack8[4], size= as1/ .pt, fontface = 'bold') +
annotate('text', x=0.6, y= -10 - 0.5, label="Birth weight", color= colorBlindBlack8[2], size= as1/ .pt, fontface = 'bold') +
  annotate('text', x=ga_fst, y=5 + 0.5, label="rs28654158", color= colorBlindBlack8[4], size= as/ .pt) +
  annotate('text', x=bw_fst, y= -10 - 0.5, label="rs11708067", color= colorBlindBlack8[2], hjust= 0, size= as/ .pt) +
  annotate('text', x= 0.6, y= 1, label= paste0('Enrichment x', round(ga_fst_pvalue, 1)), color= colorBlindBlack8[4], size= as/ .pt) +
  annotate('text', x= 0.6, y= -1, label= paste0('Enrichment x', round(bw_fst_pvalue, 1)), color= colorBlindBlack8[2], size= as/ .pt) +
  geom_density(aes(x = FST_EUR_AFR_bw, y = -..density..), fill= colorBlindBlack8[2], colour= colorBlindBlack8[2]) +
  theme_cowplot(font_size = 8) +
scale_x_continuous(expand= c(0, 0)) +
scale_y_continuous(limits= c(-11, 11), breaks= c(-10, -5, 0, 5, 10)) +
  xlab("Fst Africans - Europeans") +
ylab('Density') +
geom_segment(aes(x = ga_fst, y = 0, xend = ga_fst, yend = 5)) +
geom_segment(aes(x = bw_fst, y = 0, xend = bw_fst, yend = -10))+
geom_hline(yintercept= 0, colour= 'grey') +
theme(axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
        axis.ticks= element_line(size= 0.3))



ggsave(snakemake@output[[1]], p1, width= 63, height= 63, units= 'mm', dpi= 300)

ga_fst= d[d$v_ids== '3:123112292', 'FST_EUR_EAS']
bw_fst= d[d$v_ids== '3:123065778', 'FST_EUR_EAS']
ga_fst_pvalue= xp[xp$haplotype== 'Gestational duration' & xp$ancestries== 'FST_EUR_EAS', 'enrichment']
bw_fst_pvalue= xp[xp$haplotype== 'Birth weight' & xp$ancestries== 'FST_EUR_EAS', 'enrichment']

p1= ggplot(df1, aes(x=x) ) +
  geom_density( aes(x = FST_EUR_EAS, y = ..density..), fill= colorBlindBlack8[4], colour= 'grey') +
annotate('text', x=0.57, y= 9, label="Gestational \nduration", color= colorBlindBlack8[4], size= as1/ .pt, fontface = 'bold') +
annotate('text', x=0.57, y= -10, label="Birth weight", color= colorBlindBlack8[2], size= as1/ .pt, fontface = 'bold') +
  annotate('text', x=ga_fst, y= 10 + 0.5, label="rs28654158", color= colorBlindBlack8[4], hjust= 0, size= as/ .pt) +
  annotate('text', x=bw_fst, y= -5 - 0.5, label="rs11708067", color= colorBlindBlack8[2], size= as/ .pt) +
  annotate('text', x= 0.6, y= 1, label= paste0('Enrichment x', round(ga_fst_pvalue, 1)), color= colorBlindBlack8[4], size= as/ .pt) +
  annotate('text', x= 0.6, y= -1, label= paste0('Enrichment x', round(bw_fst_pvalue, 1)), color= colorBlindBlack8[2], size= as/ .pt) +
  geom_density( aes(x = FST_EUR_EAS_bw, y = -..density..), fill= colorBlindBlack8[2], colour= 'grey') +
scale_x_continuous(expand= c(0, 0)) +
scale_y_continuous(limits= c(-11, 11), breaks= c(-10, -5, 0, 5, 10)) +
  theme_cowplot(font_size = 8) +
  xlab("Fst East Asians - Europeans") +
ylab('Density') +
geom_segment(aes(x = ga_fst, y = 0, xend = ga_fst, yend = 10)) +
geom_segment(aes(x = bw_fst, y = 0, xend = bw_fst, yend = -5)) +
geom_hline(yintercept= 0, colour= 'grey') +
theme(axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
        axis.ticks= element_line(size= 0.3))


ggsave(snakemake@output[[2]], p1, width= 63, height= 63, units= 'mm', dpi= 300)

ga_fst= d[d$v_id== '3:123112292', 'FST_AFR_EAS']
bw_fst= d[d$v_id== '3:123065778', 'FST_AFR_EAS']
ga_fst_pvalue= xp[xp$haplotype== 'Gestational duration' & xp$ancestries== 'FST_AFR_EAS', 'enrichment']
bw_fst_pvalue= xp[xp$haplotype== 'Birth weight' & xp$ancestries== 'FST_AFR_EAS', 'enrichment']

p1= ggplot(df1, aes(x=x) ) +
geom_density( aes(x = FST_AFR_EAS, y = ..density..), fill= colorBlindBlack8[4], colour= 'grey') +
annotate('text', x=0.72, y=7, label="Gestational \nduration", color= colorBlindBlack8[4], size= as1/ .pt, fontface = 'bold') + 
annotate('text', x=0.72, y= -7, label="Birth weight", color= colorBlindBlack8[2], size= as1/ .pt, fontface = 'bold') +
annotate('text', x=ga_fst, y=5 + 0.5, label="rs28654158", color= colorBlindBlack8[4], size= as/ .pt) +
annotate('text', x=bw_fst, y= -5 - 0.5, label="rs11708067", color= colorBlindBlack8[2], hjust= 0, size= as/ .pt) +
annotate('text', x= 0.75, y= 1, label= paste0('Enrichment x', round(ga_fst_pvalue, 1)), color= colorBlindBlack8[4], size= as/ .pt) +
annotate('text', x= 0.75, y= -1, label= paste0('Enrichment x', round(bw_fst_pvalue, 1)), color= colorBlindBlack8[2], size= as/ .pt) +
geom_density( aes(x = FST_AFR_EAS_bw, y = -..density..), fill= colorBlindBlack8[2], colour= 'grey') +
scale_x_continuous(expand= c(0, 0)) +
theme_cowplot(font_size = 8) +
xlab("Fst Africans - East Asians") +
scale_y_continuous(limits= c(-11, 11), breaks= c(-10, -5, 0, 5, 10)) +
ylab('Density') +
geom_segment(aes(x = ga_fst, y = 0, xend = ga_fst, yend = 5)) +
geom_segment(aes(x = bw_fst, y = 0, xend = bw_fst, yend = -5))+
geom_hline(yintercept= 0, colour= 'grey') +
theme(axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
        axis.ticks= element_line(size= 0.3))



ggsave(snakemake@output[[3]], p1, width= 63, height= 63, units= 'mm', dpi= 300)

fwrite(df1, snakemake@output[[4]], sep= '\t')
fwrite(xp, snakemake@output[[5]], sep= '\t')
