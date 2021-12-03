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

d$ID= with(d, ifelse(REF> EFF, paste(CHR, POS, EFF, REF, sep= ':'), paste(CHR, POS, REF, EFF, sep= ':')))

outcome= ifelse(grepl('fetal', snakemake@input[[3]]), 'Fetal', 'Maternal')

x= fread(snakemake@input[[3]], select= c('ID', 'BETA', 'SE', 'pvalue'))

d= inner_join(d, x, by= 'ID')

df_h1= select(d, beta_h1, se_h1, BETA, SE)
df_h1$BETA= with(df_h1, ifelse(beta_h1<0, BETA * -1, BETA))
df_h1$beta_h1= with(df_h1, ifelse(beta_h1<0, beta_h1 * -1, beta_h1))


inputMR_m= mr_input(bx= df_h1$beta_h1, bxse= df_h1$se_h1, by= df_h1$BETA, byse= df_h1$SE)
h1= mr_allmethods(inputMR_m)$Values
names(h1)= c('method', 'estimate', 'se', 'lo95', 'up95', 'pvalue')

df_h2= select(d, beta_h2, se_h2, BETA, SE)
df_h2$BETA= with(df_h2, ifelse(beta_h2<0, BETA * -1, BETA))
df_h2$beta_h2= with(df_h2, ifelse(beta_h2<0, beta_h2 * -1, beta_h2))


inputMR_m= mr_input(bx= df_h2$beta_h2, bxse= df_h2$se_h2, by= df_h2$BETA, byse= df_h2$SE)
h2= mr_allmethods(inputMR_m)$Values
names(h2)= c('method', 'estimate', 'se', 'lo95', 'up95', 'pvalue')


df_h3= select(d, beta_h3, se_h3, BETA, SE)
df_h3$BETA= with(df_h3, ifelse(beta_h3<0, BETA * -1, BETA))
df_h3$beta_h3= with(df_h3, ifelse(beta_h3<0, beta_h3 * -1, beta_h3))

inputMR_m= mr_input(bx= df_h3$beta_h3, bxse= df_h3$se_h3, by= df_h3$BETA, byse= df_h3$SE)
h3= mr_allmethods(inputMR_m)$Values
names(h3)= c('method', 'estimate', 'se', 'lo95', 'up95', 'pvalue')

p1= ggplot(df_h1, aes(beta_h1, BETA)) +
geom_errorbarh(aes(xmin= beta_h1 - se_h1, xmax= beta_h1 + se_h1), alpha= 0.1, size= 0.3, colour= '#d9d9d9') +
geom_errorbar(aes(ymin= BETA - SE, ymax= BETA + SE), alpha= 0.1, size= 0.3, colour= '#d9d9d9') +
geom_point(color= '#737373', size= 2, alpha= 0.1, shape=21, fill= '#d9d9d9', stroke= 0.1) +
xlab('Effect of maternal transmitted\nalleles on gestational duration, days') +
ylab(paste(outcome, 'only effect\non birth weight, z-score')) +
theme_cowplot(font_size= 8) +
geom_abline(intercept= 0, slope= filter(h1, method== 'IVW') %>% pull(estimate), colour= colorBlindBlack8[2]) +
geom_abline(intercept= (filter(h1, method== '(intercept)') %>% pull(estimate))[1], slope= filter(h1, method== 'MR-Egger') %>% pull(estimate), colour= colorBlindBlack8[4]) +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.01))

p2= ggplot(df_h2, aes(beta_h2, BETA)) +
geom_errorbarh(aes(xmin= beta_h2 - se_h2, xmax= beta_h2 + se_h2), alpha= 0.1, size= 0.3, colour= '#d9d9d9') +
geom_errorbar(aes(ymin= BETA - SE, ymax= BETA + SE), alpha= 0.3, size= 0.1, colour= '#d9d9d9') +
geom_point(color= '#737373', size= 2, alpha= 0.1, fill= '#d9d9d9', shape = 21, stroke = 0.1) +
xlab('Effect of maternal non-transmitted alleles\non gestational duration, days') +
ylab(paste(outcome, 'only effect\non birth weight, z-score')) +
theme_cowplot(font_size= 8) +
geom_abline(intercept= 0, slope= filter(h2, method== 'IVW') %>% pull(estimate), colour= colorBlindBlack8[2]) +
geom_abline(intercept= (filter(h2, method== '(intercept)') %>% pull(estimate))[1], slope= filter(h2, method== 'MR-Egger') %>% pull(estimate), colour= colorBlindBlack8[4]) +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.01))

p3= ggplot(df_h3, aes(beta_h3, BETA)) +
geom_errorbarh(aes(xmin= beta_h3 - se_h3, xmax= beta_h3 + se_h3), alpha= 0.3, size= 0.1, colour= '#d9d9d9') +
geom_errorbar(aes(ymin= BETA - SE, ymax= BETA + SE), alpha= 0.3, size= 0.1, colour= '#d9d9d9') +
geom_point(color= '#737373', size= 2, alpha= 0.1, fill= '#d9d9d9', shape = 21, stroke = 0.1) +
xlab('Effect of paternal transmitted alleles\non gestational duration, days') +
ylab(paste(outcome, 'only effect\non birth weight, z-score')) +
theme_cowplot(font_size= 8) +
geom_abline(intercept= 0, slope= filter(h3, method== 'IVW') %>% pull(estimate), colour= colorBlindBlack8[2]) +
geom_abline(intercept= (filter(h3, method== '(intercept)') %>% pull(estimate))[1], slope= filter(h3, method== 'MR-Egger') %>% pull(estimate), colour= colorBlindBlack8[4]) +
geom_hline(yintercept= 0, size= 0.1) +
geom_vline(xintercept= 0, size= 0.1) +
theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks= element_blank(),
        panel.grid.major= element_line(colour= 'grey', size= 0.01))

ggsave(snakemake@output[[1]], plot= p1, width= 60, height= 60, units= 'mm', dpi= 300)
ggsave(snakemake@output[[2]], plot= p2, width= 60, height= 60, units= 'mm', dpi= 300)
ggsave(snakemake@output[[3]], plot= p3, width= 60, height= 60, units= 'mm', dpi= 300)

h1$haplotype= 'h1'
h2$haplotype= 'h2'
h3$haplotype= 'h3'

df= bind_rows(h1, h2, h3)

fwrite(d, snakemake@output[[4]], sep= '\t')
fwrite(df, snakemake@output[[5]], sep= '\t')


