library(MendelianRandomization)
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


d= fread('/mnt/work2/pol/metaGWAS/MR/raw_data/all_traits_GAraw.txt')
names(d)= c('ID', 'beta', 'se', 'pvalue', 'trait')

x= fread('/mnt/work2/pol/metaGWAS/results/meta/Maternal_GWAMA_GAraw.txt.gz', select= c('ID', 'BETA', 'SE'))

d= inner_join(d, x, by= 'ID')

funk= function(temp_df){

inputMR= mr_input(bx = temp_df$beta,   bxse= temp_df$se,by = temp_df$BETA, byse = temp_df$SE)

if (nrow(temp_df)>3) {

z= mr_allmethods(inputMR)$Values
names(z)= c('method', 'estimate', 'se', 'lo95', 'up95', 'pvalue')
z$trait= unique(temp_df$trait)

} else {
z= mr_ivw(inputMR)

z= data.frame(method= 'IVW', estimate= z$Estimate, se= z$StdError, lo95= z$CILower, up95= z$CIUpper, pvalue= z$Pvalue, trait= unique(temp_df$trait))

}
return(z)
}


mr= lapply(split(d, d$trait), funk)

mr= do.call('rbind', mr)

MR-Egger
IVW

mr= filter(mr, method== 'IVW' | method== 'MR-Egger')

mr= filter(mr, !grepl('BW', trait))

mr= arrange(mr, pvalue)

filter(mr, pvalue< 0.05/ length(unique(mr$trait)))

traits= group_by(mr, trait) %>% summarize(m1= n()) %>% filter(m1> 1) %>% pull(trait)

mr1= filter(mr, trait %in% traits)
mr2= filter(mr1, !grepl('male', trait))

mr2= arrange(mr2, desc(pvalue))

mr2$trait= factor(mr2$trait, levels= unique(mr2$trait))

ggplot(mr2, aes(trait, estimate, colour= ifelse(pvalue> 0.05 / 10 & method== 'IVW', '1', ifelse(pvalue> 0.05 / 10, '2', ifelse(method== 'IVW', '3', '4'))))) +
geom_errorbar(aes(ymin= lo95, ymax= up95), position = position_dodge(0.3), width = 1/10) +
geom_point(position = position_dodge(0.3), size= 0.8) +
scale_colour_manual(values= c('grey', 'grey',colorBlindBlack8[c(2,4)]), guide= FALSE) +
theme_cowplot(font_size= 10) +
coord_flip() +
geom_hline(yintercept= 0) +
scale_y_continuous(limits= c(-11,5)) +
ylab('Effect on gestational duration, days') +
theme(axis.title.y= element_blank())
