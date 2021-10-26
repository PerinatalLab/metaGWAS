library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
library(tidyverse)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

x= fread(snakemake@input[[1]])

x$p1= gsub('.txt.sumstats.gz', '', apply(x[, 'p1'], 1, function(x) unlist(strsplit(x, 'LDscore/'))[2]))
x$p2= gsub('.txt.sumstats.gz', '', apply(x[, 'p2'], 1, function(x) unlist(strsplit(x, 'LDSC/'))[2]))

x1= fread(snakemake@input[[2]])

x1$p1= gsub('.txt.sumstats.gz', '', apply(x1[, 'p1'], 1, function(x) unlist(strsplit(x, 'LDscore/'))[2]))
x1$p2= gsub('.txt.sumstats.gz', '', apply(x1[, 'p2'], 1, function(x) unlist(strsplit(x, 'LDSC/'))[2]))
x1$rg= -1 * x1$rg
d= rbind(x, x1)

traits= filter(d, p< 0.05/ 14, !grepl('BW', p2), !grepl('GA', p2)) %>% pull(p2)

d$trait= d$p2
d$trait= with(d, ifelse(trait== 'GAraw', 'Maternal gestational duration',
ifelse(trait== 'miscarriage', 'Miscarriage',
                ifelse(trait== 'GA_fetal', 'GA fetal effect',
                ifelse(trait== 'BW_maternal', 'Maternal BW',
                ifelse(trait== 'AFB', 'Age at first birth',
                ifelse(trait== 'AMenarche', 'Age at menarche',
                ifelse(trait== 'AMenopause', 'Age at menopause',
                ifelse(trait== 'NLB', 'Number of live births',
                ifelse(trait== 'Testosterone_fem', 'Testosterone (women)',
                ifelse(trait== 'SHBG_fem', 'SHBG (women)',
                ifelse(trait== 'SHBG_male', 'SHBG (men)',
                ifelse(trait== 'CBAT_fem', 'CBAT (women)',
                ifelse(trait== 'CBAT_male', 'CBAT (men)',
                ifelse(trait== 'Oestradiol_fem', 'Oestradiol (women)',
                ifelse(trait== 'POP', 'Pelvic Organ Prolapse',
                ifelse(trait== 'Testosterone_male', 'Testosterone (men)',
                ifelse(trait== 'leiomyoma_uterus', 'Leiomyoma uterus',
                ifelse(trait== 'BW_fetal', 'Fetal',
                ifelse(trait== 'BW_fetal_effect', 'Fetal only',
                ifelse(trait== 'Preeclampsia', 'Pre-eclampsia',
                ifelse(trait== 'BW_maternal_effect', 'Maternal only',
                ifelse(trait== 'PCOS', 'Polycistic ovary syndrome', 'Endometriosis')))))))))))))))))))))))

d= filter(d, !grepl('BW', p2), !grepl('GA', p2), !grepl('_male', p2))

traits= unique(arrange(d, p) %>% pull(trait))
d$trait= factor(d$trait, levels= rev(traits))

p1= ggplot(d, aes(rg, trait, colour= p1)) + 
geom_pointrange(aes(xmax= rg + 1.96 * se, xmin= rg - 1.96 * se), position = position_dodge(width = 0.3), fatten= 1) +
scale_colour_manual(values= colorBlindBlack8[c(8,3)], guide= FALSE) +
theme_cowplot(font_size= 9) +
scale_x_continuous(limits= c(-1, 1), breaks= seq(-1, 1, 0.5)) +
xlab('Genetic correlation') +
geom_vline(xintercept= 0, size= 0.3) +
geom_vline(xintercept= c(seq(-1, 1, 0.25)), colour= 'grey', linetype= 'dashed', alpha= 0.5, size= 0.2) +
theme(axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
        axis.ticks= element_line(size= 0.3),
        axis.title.y= element_blank())


ggsave(snakemake@output[[1]], plot= p1, width= 88, height= 120, units= 'mm', dpi= 300)
