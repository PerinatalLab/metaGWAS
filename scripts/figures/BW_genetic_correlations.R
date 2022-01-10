library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')


showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)



d= fread(snakemake@input[[1]])
d= filter(d, grepl('GAraw', p1), grepl('BW', p2))
d$p1= 'Gestational duration (maternal)'


x= fread(snakemake@input[[2]])

x= filter(x, grepl('GA_fetal', p1), grepl('BW', p2))
x$p1= 'Gestational duration (fetal)'



d= rbind(d, x)

d$p2= gsub('.txt.sumstats.gz', '', apply(d[, 'p2'], 1, function(x) unlist(strsplit(x, 'LDSC/'))[2]))
d$trait= d$p2
d$trait= with(d, ifelse(trait== 'miscarriage', 'Miscarriage',
                ifelse(trait== 'GA_fetal', 'GA fetal effect',
                ifelse(trait== 'BW_maternal', 'Maternal',
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
                ifelse(trait== 'BW_fetal_effect', 'Fetal \nonly',
                ifelse(trait== 'Preeclampsia', 'Pre-eclampsia',
                ifelse(trait== 'BW_maternal_effect', 'Maternal \nonly',
                ifelse(trait== 'PCOS', 'Polycystic ovary syndrome', 'Endometriosis'))))))))))))))))))))))



p1= ggplot(d, aes(trait, rg, colour= p1)) +
geom_pointrange(aes(ymin= rg - se * 1.96, ymax= rg + se * 1.96), position = position_dodge(0.3), width = 1/10, size= 0.4, fatten= 0.6) +
scale_colour_manual(values= colorBlindBlack8[c(8,3)], guide= FALSE) +
theme_cowplot(font_size= 8) +
scale_y_continuous(limits= c(-0.2, 0.8), breaks= seq(-0.2, 0.8, 0.2)) +
ylab('Genetic correlation') +
xlab('Effect on birth weight') +
geom_hline(yintercept= 0, size= 0.3) +
geom_hline(yintercept= c(-0.2, seq(0.2, 0.8, 0.2)), colour= 'grey', linetype= 'dashed', alpha= 0.5, size= 0.2) +
theme(axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
	axis.ticks= element_line(size= 0.3))


ggsave(snakemake@output[[1]], plot= p1, width= 60, height= 60, units= 'mm', dpi= 300)

fwrite(d, snakemake@output[[2]], sep= '\t')

