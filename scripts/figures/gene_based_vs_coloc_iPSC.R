library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
options(warn=-1)


d= fread(snakemake@input[[1]], h= T)


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)



pph= fread(snakemake@input[[1]])

geneb= fread(snakemake@input[[2]])

gene_dict= fread(snakemake@input[[3]])

names(gene_dict)= c('CHR', 'POS1', 'POS2', 'Gene', 'EnsembleID')

gene_dict$EID= with(gene_dict, unlist(lapply(strsplit(as.character(EnsembleID), ".", fixed= T), '[[', 1)))

d= inner_join(pph, gene_dict, by= c('protein'= 'EID')) %>% inner_join(., geneb, by= 'Gene')

supp_table= full_join(pph, gene_dict, by= c('protein'= 'EID')) %>% full_join(., geneb, by= 'Gene') %>% filter(Pvalue< 0.05/ nrow(geneb) | PP.H4.abf>= 0.9)


z= fread(snakemake@input[[5]], select= c('z.df1', 'z.df2', 'SNP.PP.H4', 'protein', 'snp'))
z= arrange(z, desc(SNP.PP.H4))

z= group_by(z, protein) %>% filter(row_number()==1)

d= left_join(d, z, by= c('protein'))

d= separate(d, snp, into= c('CHR', 'POS', 'REF', 'EFF'), sep= ':', remove= FALSE)

aa= fread(snakemake@input[[6]])
names(aa)= c('CHR', 'POS', 'REF', 'ALT', 'AA')
aa= filter(aa, AA!= '.')
aa= filter(aa, POS %in% d$POS)

aa$ID= with(aa, ifelse(REF> ALT, paste(CHR, POS, ALT, REF, sep= ':'), paste(CHR, POS, REF, ALT, sep= ':')))

d= left_join(d,aa[, c('ID', 'AA')], by= c('snp'= 'ID'))

d$z.df1= with(d, ifelse(d$AA== d$EFF, -1 * d$z.df1, d$z.df1))
d$z.df2= with(d, ifelse(d$AA== d$EFF, -1 * d$z.df2, d$z.df2))

d$direction= with(d, ifelse(z.df1>0 & z.df2 > 0, 'Positive', ifelse(z.df1<0 & z.df2< 0, 'Negative', 'Opposite')))
d$direction= with(d, ifelse(is.na(d$AA), 'Missing', d$direction))

#d$direction= with(d, ifelse((z.df1 * z.df2)>0, 'Same direction', 'Opposite'))

d$gene_group= with(d, ifelse(PP.H4.abf> 0.9 & Pvalue< 0.05 / nrow(geneb), 'Colocalize and gene-based significant', ifelse(Pvalue< 0.05 / nrow(geneb) & PP.H4.abf<= 0.9, 'Gene based significant',
	ifelse(PP.H4.abf> 0.9 & Pvalue> 0.05 / nrow(geneb), 'Colocalize', 'No colocalize and not significant'))))

ga= fread(snakemake@input[[4]], select= c('ID', 'BETA'))

d= inner_join(d, ga, by= c('snp'= 'ID'))

p1= ggplot(d, aes(-log10(Pvalue), PP.H4.abf, size= abs(BETA), fill= direction, alpha= (1 + PP.H4.abf) * -log10(Pvalue))) +
geom_point(shape=21, colour= 'black') +
theme_cowplot(font_size= 10) +
scale_alpha_continuous(guide= F) +
scale_size_continuous(range = c(.001, 10), guide= F) +
scale_fill_manual(values= c(colorBlindBlack8[c(1, 4, 8, 2)]), guide= F) +
geom_text_repel(data= filter(d, PP.H4.abf> 0.9 | Pvalue< 0.05 / nrow(geneb)), aes(label= Gene), max.overlaps= 20, colour= 'black', size= 6/ .pt, max.time= 10, alpha= 1) +
geom_hline(yintercept= 0.9, colour= colorBlindBlack8[8], linetype= 'dashed', size= 0.2, alpha= 0.6) +
geom_vline(xintercept= -log10(0.05/nrow(geneb)), colour= colorBlindBlack8[8], linetype= 'dashed', size= 0.2, alpha= 0.6) +
scale_y_continuous(breaks= c(seq(0, 1, 0.25), 0.9), limits= c(0, 1), expand= expansion(mult= c(0.05,0))) +
ylab('Posterior probability of colocalization') +
xlab('-log10(Gene based p-value)')

ggsave(snakemake@output[[1]], plot= p1, width= 95, height= 95, units= 'mm', dpi= 300)

d= select(d, Gene, BETA, direction, Pvalue, PP.H4.abf, Pvalue, z.df1, z.df2)

fwrite(d, snakemake@output[[2]], sep= '\t')

p1= ggplot(d, aes(-log10(Pvalue), PP.H4.abf, size= abs(BETA), fill= direction, alpha= (1 + PP.H4.abf) * -log10(Pvalue))) +
geom_point(shape=21, colour= 'black') +
theme_cowplot(font_size= 10) +
scale_alpha_continuous('Legend') +
scale_size_continuous('Legend', range = c(.001, 10)) +
scale_fill_manual('Legend', values= c(colorBlindBlack8[c(1, 4, 8, 2)])) +
geom_text_repel(data= filter(d, PP.H4.abf> 0.9 | Pvalue< 0.05 / nrow(geneb)), aes(label= Gene), max.overlaps= 20, colour= 'black', size= 6/ .pt, max.time= 10, alpha= 1) +
geom_hline(yintercept= 0.9, colour= colorBlindBlack8[8], linetype= 'dashed', size= 0.2, alpha= 0.6) +
geom_vline(xintercept= -log10(0.05/nrow(geneb)), colour= colorBlindBlack8[8], linetype= 'dashed', size= 0.2, alpha= 0.6) +
scale_y_continuous(breaks= c(seq(0, 1, 0.25), 0.9), limits= c(0, 1), expand= expansion(mult= c(0.05,0))) +
ylab('Posterior probability of colocalization') +
xlab('-log10(Gene based p-value)')

ggsave(snakemake@output[[3]], plot= p1, width= 90, height= 90, units= 'mm', dpi= 300)

fwrite(supp_table, snakemake@output[[4]], sep= '\t')
