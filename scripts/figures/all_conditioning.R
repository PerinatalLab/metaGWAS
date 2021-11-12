library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
options(warn=-1)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)


z= fread(snakemake@input[[1]])
z$chr= as.numeric(gsub('chr', '', z$chr))
z$chr= as.character(z$chr)
z$locus= 1:nrow(z)


funk= function(infile){
d= fread(infile)
names(d)[1:11]= names(d)[2:12]
d=d[, 1:11]

d= filter(d, p<5e-6)

d$bC= ifelse(d$b< 0, -1 * d$bC, d$bC)
d$b= ifelse(d$b< 0, -1 * d$b, d$b)

d= separate(d, SNP, into= c('chr', 'POS', 'REF', 'EFF'), sep= ':')

d$POS= as.numeric(d$POS)
d$chr= as.character(d$chr)
d$GWAS= ifelse(grepl('BW_maternal_effect_GA', infile), 'BW_maternal_GA', ifelse(grepl('BW_fetal_effect_GA', infile), 'BW_fetal_GA', ifelse(grepl('GA_BW_maternal', infile), 'GA_BW_maternal', 'GA_BW_fetal')))
d= inner_join(d, z, on= 'chr') 
d= d %>% filter(POS>= start, POS< stop)

d= group_by(d, locus) %>% arrange(p) %>% filter(row_number()== 1)

return(d)

}

df_list= lapply(snakemake@input[grepl('BW', snakemake@input)], funk)

d= do.call('rbind', df_list)

d$beta_dif= with(d, (bC - b) / b)




p1= ggplot(d, aes(GWAS, beta_dif, fill= GWAS)) + 
geom_boxplot() +
theme_cowplot(font_size= 10) +
scale_fill_manual(values= c(colorBlindBlack8[c(4, 2, 8, 7)], 'grey'), guide= F) +
geom_jitter(size=0.4, alpha=0.4) +
xlab('\nP-value on birth weight (maternal only effect)') + 
ylab('Relative difference in effect size after \nconditioning birth weight on gestational duration') +
geom_hline(yintercept= 0, linetype= 'dashed', alpha= 0.6, colour= 'grey')

ggsave(snakemake@output[[1]], plot= p1, width= 90, height= 90, units= 'mm', dpi= 300)


fwrite(d, snakemake@output[[2]], sep= '\t')


