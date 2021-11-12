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

ge= data.frame(CHR= c(5, 3, 1, 23, 3, 23), pos_ge= c(157895049, 127881613, 22470407, 115164770, 123068359, 131300571))

df= arrange(d, pvalue)

df= df[!duplicated(df[, c('CHR', 'POS')]), ]

dg= fread(snakemake@input[[2]])
dg$GENE= dg$nearestGene


don <- df %>%
    group_by(CHR)      %>%
    summarise(chr_len= max(POS)) %>%
    mutate(tot= cumsum(as.numeric(chr_len))-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(df, ., by= 'CHR') %>%
    arrange(CHR, POS) %>% # Add a cumulative position of each SNP
    mutate( BPcum=POS+tot) %>%
         ungroup()

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  names(axisdf)= c('CHR', 'center')
HC= -log10(5*10**-8)
dg= dg %>% ungroup() %>% select(ID, GENE, CHR, POS, BETA)
don= left_join(don, select(dg, ID, GENE), by= 'ID')
names(dg)= c('ID', 'GENE', 'CHR', 'POS_new', 'BETA')

lims= 250000

don= data.frame(don)
dg= data.frame(dg)

don$disc= ifelse(don$pvalue> 5*10**-8, 0, 2)

for (i in rownames(dg)) {
don= mutate(don, disc= ifelse(CHR== as.integer(dg[i, 'CHR']) & POS>= as.integer(dg[i, 'POS_new']) - lims & POS<= as.integer(dg[i, 'POS_new']) + lims, 2, disc))
}

for (i in rownames(ge)) {
don= mutate(don, disc= ifelse(CHR== as.integer(ge[i, 'CHR']) & POS>= as.integer(ge[i, 'pos_ge']) - lims & POS<= as.integer(ge[i, 'pos_ge']) + lims, 1, disc))
}

don= don[order(don$disc, decreasing= F, na.last= T), ]
don$disc= factor(don$disc, levels=c(0, 1, 2), labels=c('Not significant', 'Previous discovery', 'New discovery'))

cols <- c('Not significant'= 'grey', 'Previous discovery'= colorBlindBlack8[4], 'New discovery'= colorBlindBlack8[2])

don$GENE= ifelse(!is.na(don$GENE), don$nearestGene, don$GENE)

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)


 p1= ggplot(data= don, aes(x=-log10(pvalue), y= BPcum, colour= disc)) +
  geom_point(size= 0.07) +   # Show all points
  theme_cowplot(font_size= 9) +
  #scale_alpha_manual(values= rep(c(1/10, 1/2), 23)) +
  scale_colour_manual(values= cols, guide= F) +
  scale_y_reverse(label = ifelse(axisdf$CHR== 23, 'X', axisdf$CHR), breaks= axisdf$center, expand= c(0.03, 0.03)) + # custom y axis
  scale_x_continuous(expand= c(0, 0), limits= c(0, 43), breaks= seq(0, 40, 10), sec.axis = sec_axis(~ ., name = derive())) +
  xlab('-log10(pvalue)') +
  ylab('Chromosome') +
  geom_vline(xintercept= 0, size= 0.5, colour= 'black') +
  geom_vline(xintercept= HC, size= 0.5, linetype= 2, colour= '#878787') +
  coord_cartesian(clip = "off") +
  geom_text_repel(data= filter(don, GENE!= ''), aes(y= BPcum, x= -log10(pvalue), label= GENE),
                  size= 6/ .pt,
                  force_pull= 0, # do not pull toward data points
                  force= 0.1,
                  nudge_x      = 43 - ((-log10(filter(don, GENE!= '')$pvalue))),
                  direction    = "y",
                  angle        = 00,
                  hjust        = 0,
                  vjust=  0.5,
		  box.padding= 0.1,
                  segment.size = 0.1,
                  segment.square= TRUE,
                  segment.inflect= FALSE,
                  segment.colour= colorBlindBlack8[8],
                  colour= 'black',
                  segment.linetype = 4,
                  xlim = c(-Inf, 50),
                  ylim = c(-Inf, Inf)) +
  theme(legend.position= 'none',
	plot.margin = unit(c(t= 0, r=1, b= 0, l=0), 'cm'),
        text= element_text(family="arial", size= 9),
	axis.line= element_line(size= 0.1) )  

save_plot(snakemake@output[[1]], plot= p1, base_height= 160, base_width= 90, units= 'mm', dpi= 300)

