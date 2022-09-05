library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
options(warn=-1)


d= fread(snakemake@input[[1]], h= T, select= c('ID', 'CHR', 'POS', 'pvalue', 'nearestGene'))
d$pheno= 'GAraw'
x= fread(snakemake@input[[3]], h= T, select= c('ID', 'CHR', 'POS', 'pvalue', 'nearestGene'))
x$pheno= 'allPTD'

d= rbind(d, x)

rm(x)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ge= data.frame(CHR= c(5, 3, 1, 23, 3, 23), pos_ge= c(157895049, 127881613, 22470407, 115164770, 123068359, 131300571))

df= arrange(d, pvalue)


dg= fread(snakemake@input[[2]])
dg$GENE= dg$nearestGene

ptd= fread(snakemake@input[[4]])
ptd$GENE= ptd$nearestGene

don <- df %>%
    group_by(CHR)      %>%
    summarise(chr_len= max(POS)) %>%
    mutate(tot= cumsum(as.numeric(chr_len))-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(df, ., by= 'CHR') %>%
    arrange(CHR, POS) %>% # Add a cumulative position of each SNP
    mutate(BPcum=POS+tot) %>%
         ungroup()

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  names(axisdf)= c('CHR', 'center')
HC= -log10(5*10**-8)
dg= dg %>% ungroup() %>% select(ID, GENE, CHR, POS, BETA)
ptd= ptd %>% ungroup %>% select(ID, GENE, CHR, POS, BETA)

don$disc= ifelse(don$pvalue> 5*10**-8, 0, 2)

don1= filter(don, pheno== 'GAraw') %>% left_join(., select(dg, ID, GENE), by= 'ID')
don2= filter(don, pheno!= 'GAraw') %>% left_join(., select(ptd, ID, GENE), by= 'ID')
names(dg)= c('ID', 'GENE', 'CHR', 'POS_new', 'BETA')
names(ptd)= c('ID', 'GENE', 'CHR', 'POS_new', 'BETA')

lims= 250000

don= data.frame(don)
dg= data.frame(dg)
ptd= data.frame(ptd)


for (i in rownames(dg)) {
don1= mutate(don1, disc= ifelse(CHR== as.integer(dg[i, 'CHR']) & POS>= as.integer(dg[i, 'POS_new']) - lims & POS<= as.integer(dg[i, 'POS_new']) + lims, 2, disc)) 
}

for (i in rownames(ptd)) {
don2= mutate(don2, disc= ifelse(CHR== as.integer(ptd[i, 'CHR']) & POS>= as.integer(ptd[i, 'POS_new']) - lims & POS<= as.integer(ptd[i, 'POS_new']) + lims, 2, disc))

}

don= rbind(don1, don2)
rm(don1) ; rm(don2)

for (i in rownames(ge)) {
don= mutate(don, disc= ifelse(CHR== as.integer(ge[i, 'CHR']) & POS>= as.integer(ge[i, 'pos_ge']) - lims & POS<= as.integer(ge[i, 'pos_ge']) + lims, 1, disc))
}

don= don[order(don$disc, decreasing= F, na.last= T), ]
don$disc= factor(don$disc, levels=c(0, 1, 2), labels=c('Not significant', 'Previous discovery', 'New discovery'))

cols <- c('Not significant'= 'grey', 'Previous discovery'= colorBlindBlack8[3], 'New discovery'= colorBlindBlack8[8])

don$GENE= ifelse(!is.na(don$GENE), don$nearestGene, don$GENE)

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)


don$GENE= with(don, ifelse(GENE== 'CDC42', 'CDC42/ WNT4', ifelse(GENE== 'HIVEP3', 'HIVEP3/ EDN2', ifelse(GENE== 'TET3', 'TET3/ DGUOK-AS1', ifelse(GENE== 'TCEA2', 'TCEA2/ OPRL1', GENE)))))

don$logpval= with(don, ifelse(pheno== 'allPTD', log10(pvalue), -log10(pvalue)))

p1= ggplot(data= don, aes(x= BPcum, y= logpval, colour= disc)) +
  geom_point(size= 0.07) +   # Show all points
  theme_cowplot(font_size= 9) +
  scale_colour_manual(values= cols, guide= F) +
  scale_x_continuous(label = c(1:19, '', 21,'', 'X'), breaks= axisdf$center, expand= c(0.03, 0.03)) + # label = ifelse(axisdf$CHR== 23, 'X', axisdf$CHR)
  scale_y_continuous(expand= c(0, 0), limits= c(min(don$logpval) - 2, max(don$logpval) + 2), breaks= seq(-30, 45, 10), labels= c(abs(seq(-30, 45, 10)))) + # , sec.axis = sec_axis(~ ., name = derive())) +
  ylab('-log10(pvalue)') +
  xlab('Chromosome') +
  geom_hline(yintercept= 0,, size= 0.25, colour= 'black') +
  geom_hline(yintercept= c(HC, -HC), size= 0.2, linetype= 2, colour= '#878787') +
  coord_cartesian(clip = "off") +
  geom_text_repel(data= filter(don, GENE!= ''), aes(x= BPcum, y= logpval, label= GENE),
                  size= 6/ .pt,
                  force_pull= 0, # do not pull toward data points
                  force= 0.1,
                  nudge_y      =  ifelse(filter(don, GENE!= '') %>% pull(logpval)>0, 1, -1), #43 - ((-log10(filter(don, GENE!= '')$pvalue))),
                  direction    = "both",
                  hjust        = 0,
                  vjust=  0.5,
		  box.padding= 0.1,
		  angle= 0,
                  segment.size = 0.1,
                  segment.square= TRUE,
                  segment.inflect= FALSE,
                  segment.colour= colorBlindBlack8[8],
                  colour= ifelse(filter(don, GENE!= '') %>% pull(disc)== 'New discovery', colorBlindBlack8[8], colorBlindBlack8[3]),
                  segment.linetype = 4,
                  ylim = c(-Inf, 50),
                  xlim = c(-Inf, Inf)) +
  theme(legend.position= 'none',
	plot.margin = unit(c(t= 0, r=0, b= 0, l=0), 'cm'),
        text= element_text(family="arial", size= 9),
	axis.line= element_line(size= 0.1)) 

save_plot(snakemake@output[[1]], plot= p1, base_height= 90, base_width= 180, units= 'mm', dpi= 300)

