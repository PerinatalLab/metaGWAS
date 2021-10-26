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

d= filter(d, grepl('GAraw', p1), !grepl('BW', p2), !grepl('male', p2))
#d$p1= 'Gestational duration (maternal)'
d$p1= 'Gestational duration - Maternal'
x= fread(snakemake@input[[2]])

bw= filter(x, grepl('BW_maternal_effect', p1), !grepl('BW', p2), !grepl('male', p2), !grepl('GA_fetal', p2))
x= filter(x, grepl('GA_fetal', p1), !grepl('BW', p2), !grepl('male', p2))
#x$p1= 'Gestational duration (fetal)'
x$p1= 'Gestational duration - Fetal'

bw$p1= 'Birth weight - Maternal'

bwc= fread(snakemake@input[[3]])
bwc$p1= 'Birth weight conditioned - Maternal'
bwc= filter(bwc, !grepl('BW', p2), !grepl('male', p2), !grepl('GA_fetal', p2))

d= do.call('rbind', list(d, x, bw,bwc))

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
                ifelse(trait== 'BW_fetal_effect', 'Fetal only',
                ifelse(trait== 'Preeclampsia', 'Pre-eclampsia',
                ifelse(trait== 'BW_maternal_effect', 'Maternal only',
                ifelse(trait== 'PCOS', 'Polycistic ovary syndrome', 'Endometriosis'))))))))))))))))))))))



d= filter(d, trait!= 'GA fetal effect')


traits_list= filter(d, p1== 'Gestational duration - Maternal') %>% arrange(rg) %>% pull(trait)
d$trait= factor(d$trait, levels= unique(traits_list))

d$sig= ifelse(d$p< 0.05/ (nrow(d)/2), '**', ifelse(d$p< 0.05, '*', ''))

d$p1= factor(d$p1, levels= rev(c('Gestational duration - Maternal', 'Birth weight - Maternal', 'Birth weight conditioned - Maternal', 'Gestational duration - Fetal')))

p1= ggplot(d, aes(trait, p1, fill= rg)) +
geom_tile(colour = "white", size= 1) +
theme_cowplot(font_size= 9) +
scale_fill_gradient2(low= colorBlindBlack8[2], high= colorBlindBlack8[4], mid= 'white', guide= F) +
theme(axis.text.x = element_text(angle = 45, hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
	axis.ticks.y= element_blank()) +
scale_x_discrete(position = "top") +
geom_text(data= d, aes(trait, p1, label= sig), size= 6/ .pt) +
theme(  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour= 'black', fill= NA, size=1),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        axis.line= element_blank()) +
coord_equal()



ggsave(snakemake@output[[1]], plot= p1, width= 127, height= 85, units= 'mm', dpi= 300)

fwrite(d, snakemake@output[[2]], sep= '\t')


p1= ggplot(d, aes(trait, p1, fill= rg)) +
geom_tile(colour = "white", size= 1) +
theme_cowplot(font_size= 10) +
scale_fill_gradient2(low= colorBlindBlack8[2], high= colorBlindBlack8[4], mid= 'white') +
theme(axis.text.x= element_text(angle = 45, hjust = 0),
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        plot.margin= margin(t= 7, r= 32, l= 1, b= 7, unit= 'pt')) +
scale_x_discrete(position = "top") +
geom_text(data= d, aes(trait, p1, label= sig), size= 6/ .pt)

ggsave(snakemake@output[[3]], plot= p1, width= 127, height= 45, units= 'mm', dpi= 300)

