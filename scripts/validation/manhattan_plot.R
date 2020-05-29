library(data.table)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(cowplot)
library('ggrepel')

colors_2= c('#9C02A7', '#FFBD01')

d= fread(snakemake@input[[1]])

gene= fread(snakemake@input[[3]])

df= fread(snakemake@input[[2]])
df= separate(df, variant, into= c('chr', 'pos', 'ref', 'eff'), ':')

df$chr= gsub('X', '23', df$chr)
df$chr= as.numeric(df$chr)
df$pos= as.numeric(df$pos)

df= filter(df, low_confidence_variant== F)

d= full_join(d, df, by= c('chr', 'pos'))

df= d

df= df %>% filter(!is.na(chr))

don <- df %>%
    group_by(chr)      %>%
    summarise(chr_len= max(pos)) %>%
    mutate(tot= cumsum(chr_len)-chr_len) %>% # Calculate cumulative position of each chromosome
    select(-chr_len) %>%
    left_join(df, ., by= 'chr') %>%
    arrange(chr, pos) %>% # Add a cumulative position of each SNP
    mutate( BPcum=pos+tot)

axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
names(axisdf)= c('chr', 'center')



red_don= filter(don, metaP< 10**-4)

gene= inner_join(select(red_don, chr, BPcum, pos, metaP), gene, by= 'chr')

gene= filter(gene, pos> start- 10000, pos< end + 10000)

gene= arrange(gene, chr, pos) %>% group_by(chr) %>% mutate(d=pos-lag(pos, default=-Inf), clumpid=cumsum(d>250000)) %>% group_by(chr, clumpid) %>% filter(rank(metaP, ties.method = "random")==1)

don$replicated= ifelse(!is.na(don$metaP) & don$pval<10**-4, 'rep', 'not')

don1= filter(don, pval> 0.05, is.na(metaP))
don= filter(don, pval<= 0.05, !is.na(metaP))
don1$pval= round(don1$pval, 3)
don1= distinct(don1, chr, pval, .keep_all=T)

don= rbind(don, don1)

p1= ggplot(don) +
geom_point(data= don, aes(x=BPcum, y= -log10(metaP), color= factor(replicated), alpha=  factor(chr)), size=0.3, shape= 21) +   # Show all points +
geom_point(data= don, aes(x=BPcum, y= log10(pval), color= factor(replicated), alpha= factor(chr)),size=0.3, shape= 21) +
theme_cowplot(12, font_size= 12) + #theme_minimal_hgrid(12, rel_small = -1) + 
scale_color_manual(values= c('rep'= '#9C02A7', 'not'= 'grey')) +
scale_alpha_manual(values= rep(c(1, 0.5), 12)) +
scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, expand=c(0,0) ) + # custom X axis
theme( legend.position="none") +
xlab('Chromosome') +
ylab('-/ log10(p-value)') +
geom_hline(yintercept= 0, size= 0.5, colour= 'black') +
geom_hline(yintercept= -log10(5*10**-8), size= 0.5, linetype= 2, colour= '#878787') +
geom_hline(yintercept= log10(10**-5), size= 0.5, linetype= 2, colour= '#878787') +
#annotate(geom="text", x= Inf, y= HC - 0.5, label= 'bold("High confidence")', color="black", vjust= 1, hjust= 1, parse= TRUE, size= 4)
geom_text_repel(data= gene, aes(x= BPcum, y= -log10(metaP), label= gene), size= 3,  hjust = 0.5, force= 1, vjust= 1, colour= 'black')


save_plot(file= snakemake@output[[1]], plot= p1, base_width=297, base_height=210, units="mm", device= cairo_ps)
