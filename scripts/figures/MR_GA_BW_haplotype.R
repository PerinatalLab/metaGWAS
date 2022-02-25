library(MendelianRandomization)
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

d$rsid= with(d, ifelse(rsid== 'chrX:116013571', 'rs5991030', ifelse(rsid== 'chrX:132178061', 'rs5930554', rsid)))

d$effect_origin= with(d, ifelse(class_name== 'MF OD' | class_name== 'MF SD', 'Maternal and fetal', ifelse(class_name== 'Fetal MatT' | class_name== 'Fetal', 'Fetal', 'Maternal')))

#d= filter(d, MarkerName!= '6:32595083:G:T')

#top= fread(snakemake@input[[2]])
#ids= pull(top, ID)
#ids= c('3:156697097:A:G', '5:158058432:G:T', ids)

x= fread(snakemake@input[[2]], select= c('ID', 'RSID'))

#x= filter(x, ID %in% ids)

d= inner_join(d, x, by= c('rsid'= 'RSID'))

d= separate(d, ID, into= c('CHR', 'POS', 'REF', 'EFF'), sep= ':')
d$beta_MT= with(d, ifelse(REF > EFF, -1 * beta_MT, beta_MT))
d$beta_MNT= with(d, ifelse(REF > EFF, -1 * beta_MNT, beta_MNT))
d$beta_PT= with(d, ifelse(REF > EFF, -1 * beta_PT, beta_PT))

d$ID= with(d, ifelse(REF> EFF, paste(CHR, POS, EFF, REF, sep= ':'), paste(CHR, POS, REF, EFF, sep= ':')))

outcome= ifelse(grepl('fetal', snakemake@input[[3]]), 'Fetal', 'Maternal')

x= fread(snakemake@input[[3]], select= c('ID', 'BETA', 'SE', 'pvalue'))

d= inner_join(d, x, by= 'ID')

df_MT= select(d, beta_MT, se_MT, BETA, SE, effect_origin)
df_MT$BETA= with(df_MT, ifelse(beta_MT<0, BETA * -1, BETA))
df_MT$beta_MT= with(df_MT, ifelse(beta_MT<0, beta_MT * -1, beta_MT))


inputMR_m= mr_input(bx= df_MT$beta_MT, bxse= df_MT$se_MT, by= df_MT$BETA, byse= df_MT$SE)
MT= mr_allmethods(inputMR_m)$Values
names(MT)= c('method', 'estimate', 'se', 'lo95', 'up95', 'pvalue')

df_MNT= select(d, beta_MNT, se_MNT, BETA, SE, effect_origin)
df_MNT$BETA= with(df_MNT, ifelse(beta_MNT<0, BETA * -1, BETA))
df_MNT$beta_MNT= with(df_MNT, ifelse(beta_MNT<0, beta_MNT * -1, beta_MNT))


inputMR_m= mr_input(bx= df_MNT$beta_MNT, bxse= df_MNT$se_MNT, by= df_MNT$BETA, byse= df_MNT$SE)
MNT= mr_allmethods(inputMR_m)$Values
names(MNT)= c('method', 'estimate', 'se', 'lo95', 'up95', 'pvalue')


df_PT= select(d, beta_PT, se_PT, BETA, SE, effect_origin)
print(nrow(df_PT))
df_PT$BETA= with(df_PT, ifelse(beta_PT<0, BETA * -1, BETA))
df_PT$beta_PT= with(df_PT, ifelse(beta_PT<0, beta_PT * -1, beta_PT))

inputMR_m= mr_input(bx= df_PT$beta_PT, bxse= df_PT$se_PT, by= df_PT$BETA, byse= df_PT$SE)
PT= mr_allmethods(inputMR_m)$Values
names(PT)= c('method', 'estimate', 'se', 'lo95', 'up95', 'pvalue')

p1= ggplot(df_MT, aes(beta_MT, BETA, colour= effect_origin, fill= effect_origin)) +
geom_errorbarh(aes(xmin= beta_MT - se_MT, xmax= beta_MT + se_MT, colour= effect_origin, fill= effect_origin), size= 0.1, alpha= 0.7) +
geom_errorbar(aes(ymin= BETA - SE, ymax= BETA + SE, colour= effect_origin, fill= effect_origin),size= 0.1, alpha= 0.7) +
geom_point(size= 2, shape=21, stroke= 0.1, alpha= 0.7) +
scale_colour_manual(values= colorBlindBlack8[c(4, 2, 1)], guide= 'none') +
scale_fill_manual(values= colorBlindBlack8[c(4, 2, 1)], guide= 'none') +
xlab('Effect of maternal transmitted\nalleles on gestational duration, days') +
ylab(paste(outcome, 'only effect\non birth weight, z-score')) +
theme_cowplot(font_size= 8) +
geom_abline(intercept= 0, slope= filter(MT, method== 'IVW') %>% pull(estimate), colour= '#d9d9d9') +
geom_abline(intercept= (filter(MT, method== '(intercept)') %>% pull(estimate))[1], slope= filter(MT, method== 'MR-Egger') %>% pull(estimate), colour= '#d9d9d9', linetype= 'dashed') +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.05))

p2= ggplot(df_MNT, aes(beta_MNT, BETA, colour= effect_origin, fill= effect_origin)) +
geom_errorbarh(aes(xmin= beta_MNT - se_MNT, xmax= beta_MNT + se_MNT,colour= effect_origin, fill= effect_origin), size= 0.1) +
geom_errorbar(aes(ymin= BETA - SE, ymax= BETA + SE,colour= effect_origin, fill= effect_origin),size= 0.1) +
geom_point(size= 2, shape= 21, stroke= 0.1) +
scale_colour_manual(values= alpha(colorBlindBlack8[c(4, 2, 1)], 0.7), guide= 'none') +
scale_fill_manual(values= alpha(colorBlindBlack8[c(4, 2, 1)], 0.7), guide= 'none') +
xlab('Effect of maternal non-transmitted alleles\non gestational duration, days') +
ylab(paste(outcome, 'only effect\non birth weight, z-score')) +
theme_cowplot(font_size= 8) +
geom_abline(intercept= 0, slope= filter(MNT, method== 'IVW') %>% pull(estimate), colour= '#d9d9d9') +
geom_abline(intercept= (filter(MNT, method== '(intercept)') %>% pull(estimate))[1], slope= filter(MNT, method== 'MR-Egger') %>% pull(estimate), colour= '#d9d9d9', linetype= 'dashed') +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.05))

p3= ggplot(df_PT, aes(beta_PT, BETA, colour= effect_origin, fill= effect_origin)) +
geom_errorbarh(aes(xmin= beta_PT - se_PT, xmax= beta_PT + se_PT, colour= effect_origin, fill= effect_origin), size= 0.1) +
geom_errorbar(aes(ymin= BETA - SE, ymax= BETA + SE, colour= effect_origin, fill= effect_origin), alpha= 0.5, size= 0.1) +
geom_point(size= 2, shape= 21, stroke = 0.1) +
scale_colour_manual(values= alpha(colorBlindBlack8[c(4, 2, 1)], 0.7), guide= 'none') +
scale_fill_manual(values= alpha(colorBlindBlack8[c(4, 2, 1)], 0.7), guide= 'none') +
xlab('Effect of paternal transmitted alleles\non gestational duration, days') +
ylab(paste(outcome, 'only effect\non birth weight, z-score')) +
theme_cowplot(font_size= 8) +
geom_abline(intercept= 0, slope= filter(PT, method== 'IVW') %>% pull(estimate), colour= '#d9d9d9') +
geom_abline(intercept= (filter(PT, method== '(intercept)') %>% pull(estimate))[1], slope= filter(PT, method== 'MR-Egger') %>% pull(estimate), colour= '#d9d9d9', linetype= 'dashed') +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.05))

ggsave(snakemake@output[[1]], plot= p1, width= 70, height= 70, units= 'mm', dpi= 300)
ggsave(snakemake@output[[2]], plot= p2, width= 70, height= 70, units= 'mm', dpi= 300)
ggsave(snakemake@output[[3]], plot= p3, width= 70, height= 70, units= 'mm', dpi= 300)

MT$haplotype= 'MT'
MNT$haplotype= 'MNT'
PT$haplotype= 'PT'

df= bind_rows(MT, MNT, PT)

fwrite(d, snakemake@output[[4]], sep= '\t')
fwrite(df, snakemake@output[[5]], sep= '\t')


