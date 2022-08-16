library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
options(warn=-1)

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)


colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

d= fread(snakemake@input[[1]])

z= fread(snakemake@input[[3]])

df= fread(snakemake@input[[2]], select= (c('MarkerName', 'Effect', 'StdErr', 'HetISq', 'HetPVal', 'TOTALSAMPLESIZE', 'P-value', 'Allele1', 'Allele2')))
names(df)= c('SNP', 'BETA', 'SE', 'HetISq', 'HetPval', 'N', 'pvalue', 'A1', 'A2')
df= filter(df, SNP %in% d$SNP)

df= separate(df, SNP, into= c('CHR', 'POS', 'Ax1', 'Ax2', 'ID'), sep= ':', remove= F)
df$BETA= ifelse(df$A2 > df$A1, -1 * df$BETA, df$BETA)
df$CHR= ifelse(df$CHR== 'X','23', df$CHR)
df$CHR= as.integer(df$CHR)
df$POS= as.integer(df$POS)
df= select(df, -c(A1, A2, ID, Ax1, Ax2))

df$cohort= 'Meta-analysis'
d= bind_rows(d, df)

z$CHR= ifelse(z$CHR== 'X','23', z$CHR)
z$CHR= as.integer(z$CHR)

d= inner_join(d, z, by= 'CHR') %>% filter(POS> pos1, POS< pos2)

d$locus= paste0('Chr ', d$CHR,': ', d$nearestGene)

d= filter(d, !(cohort %in% c('PGPII', 'PGPIII', 'BIB', 'DNBCPTD', 'STORK', 'STORKGROR')))

d$cohort= paste0(d$cohort, ' (n= ', d$N, ')')

temp_df= d[d$nearestGene== 'EEFSEC', ]

temp_df= temp_df[order(temp_df$N, decreasing= T), ]

p1= ggplot(temp_df, aes(x=factor(cohort, level = factor(cohort)), y=BETA, ymin= BETA - 1.96 * SE, ymax= BETA + 1.96 * SE, colour= !is.na(HetISq), shape= !is.na(HetISq)), alpha= 0.5) +
 geom_pointrange(size= 1) +
scale_shape_manual(values= c(15, 18), guide= F) +
 geom_hline(yintercept = 0, linetype=2) +
scale_y_continuous(sec.axis = dup_axis()) +
 coord_flip() +
scale_colour_manual(values= c(colorBlindBlack8[3], colorBlindBlack8[4]), guide= F) +
theme_cowplot() +
 xlab('') +
    ylab('Beta [95% CI]') +
geom_vline(xintercept= 0, linetype= "dotted", colour= 'grey') 

ggsave(snakemake@output[[1]], plot= p1, width= 120, height= 90, units= 'mm', dpi= 300)

