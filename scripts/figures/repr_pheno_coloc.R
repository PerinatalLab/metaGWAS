library(data.table)
library(dplyr)
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library('showtext')


colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

d= fread(snakemake@input[[1]])


d= filter(d, !grepl('BW', trait), !grepl('GA_fetal', trait), !grepl('male', trait))

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
                ifelse(trait== 'PCOS', 'Polycystic ovary syndrome', 'Endometriosis'))))))))))))))))))))))

pregnancy= c('Miscarriage', 'Pre-eclampsia')
uterus= c('Leiomyoma uterus', 'Pelvic Organ Prolapse', 'Endometriosis', 'Polycystic ovary syndrome')
fitness= c('Age at first birth', 'Number of live births')
hormonal= c('Age at menarche', 'Age at menopause', 'Testosterone (women)', 'SHBG (women)', 'CBAT (women)', 'Oestradiol (women)')

d$cluster= with(d, ifelse(trait %in% pregnancy, 'Pregnancy', ifelse(trait %in% uterus, 'Reproductive organs', ifelse(trait %in% fitness, 'Fitness', 'Sex-hormone related'))))

d$colour= with(d, ifelse(cluster== 'Pregnancy', colorBlindBlack8[3], ifelse(cluster== 'Reproductive organs', colorBlindBlack8[5], ifelse(cluster== 'Fitness', colorBlindBlack8[7], colorBlindBlack8[8]))))

d$GENE= apply(d[, 'locus'], 1, function(x) unlist(strsplit(x, '_'))[2])

d$GENE= with(d, ifelse(GENE== 'CDC42', 'CDC42/ WNT4', ifelse(GENE== 'HIVEP3', 'HIVEP3/ EDN2', ifelse(GENE== 'TET3', 'TET3/ DGUOK-AS1', ifelse(GENE== 'TCEA2', 'TCEA2/ OPRL1', GENE)))))

d$sig= ifelse(d$PP.H4.abf>0.5, '*', '')

d= arrange(d, cluster)

d$trait= factor(d$trait, levels= unique(d$trait))
traits= unique(d$trait)
colors <- filter(d, !duplicated(trait)) %>% arrange(trait) %>% pull(colour)

d$PP= ifelse(d$PP.H4.abf> d$PP.H3.abf, d$PP.H4.abf, -d$PP.H3.abf - d$PP.H4.abf)
d$PP2= ifelse(d$PP.H4.abf> d$PP.H3.abf, d$PP.H4.abf, d$PP.H3.abf)
p1= ggplot(d, aes(trait, GENE, value= PP, fill= PP, colour= PP, size= PP2, stroke= 1-  PP)) + 
theme_cowplot(font_size= 9) +
geom_point(shape= 15) + 
scale_fill_gradient2(low= colorBlindBlack8[4], mid= 'white', high= colorBlindBlack8[2], guide= F) +
scale_colour_gradient2(low= colorBlindBlack8[4], mid= 'white', high= colorBlindBlack8[2], guide= F) +
scale_size_continuous(range= c(1, 2.5), guide= F) +
scale_x_discrete(position= 'top') +
theme(axis.ticks= element_blank(),
	axis.title= element_blank(),
	axis.text.x= element_blank())  +
geom_vline(xintercept= 1:(length(unique(d$trait))-1) + 0.5, size= 0.4, colour= 'grey') +
geom_hline(yintercept= 1:(length(unique(d$GENE))-1) + 0.5, size= 0.4, colour= 'grey') +
geom_vline(xintercept= cumsum(c(length(fitness) , length(pregnancy) , length(uterus) )) +0.5, size= 0.8) +
theme(	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	panel.background = element_blank(),
	panel.border = element_rect(colour= 'black', fill= NA, size=1),
	plot.margin = unit(c(0, 0.1, 0.1, 0), "cm"),
	axis.line= element_blank())

t_count_locus= group_by(d, trait) %>% summarize(PP= sum(as.numeric(PP.H4.abf> 0.8)), PP_locus= sum(as.numeric(PP.H4.abf + PP.H3.abf>0.8)))
t_count_locus$PP= t_count_locus$PP_locus - t_count_locus$PP
t_count_locus$supp= 'Locus-level'

t_count= group_by(d, trait) %>% summarize(PP= sum(as.numeric(PP.H4.abf> 0.8)))
t_count$supp= 'Coloc'

t_count= bind_rows(t_count, t_count_locus)

t_count$trait= factor(t_count$trait, levels= unique(d$trait))
t_count$supp= factor(t_count$supp, levels= c('Locus-level','Coloc'))

p2= ggplot(t_count, aes(trait, -PP, fill= supp)) +
theme_cowplot(font_size= 8) +
geom_col(alpha= 0.7) +
geom_hline(yintercept= 0) +
scale_fill_manual(values= c(colorBlindBlack8[4], colorBlindBlack8[2]), guide= F) +
theme(	axis.line= element_blank(),
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	panel.background = element_blank(),
	panel.border = element_rect(colour= 'black', fill= NA, size=1),
	axis.text.x= element_blank(),
	axis.ticks.x= element_blank(),
	axis.title= element_blank(),
	plot.margin = unit(c(0, 0, 0, 0.1), "cm")) +
scale_y_continuous(limits= c(-10, 0), expand= c(0,0), labels= seq(0, 10, 2), breaks= seq(0, -10, -2)) +
geom_vline(xintercept= cumsum(c(length(fitness) , length(pregnancy) , length(uterus) )) +0.5, size= 0.8)  +
geom_hline(yintercept= c(-4, -8), size= 0.3, linetype= 'dashed', colour= 'grey')

l_count_locus= group_by(d, GENE) %>% summarize(PP= sum(as.numeric(PP.H4.abf> 0.8)), PP_locus= sum(as.numeric(PP.H4.abf + PP.H3.abf>0.8)))
l_count_locus$PP= l_count_locus$PP_locus - l_count_locus$PP
l_count_locus$supp= 'Locus-level'

l_count= group_by(d, GENE) %>% summarize(PP= sum(as.numeric(PP.H4.abf> 0.8)))
l_count$supp= 'Coloc'

l_count= bind_rows(l_count, l_count_locus)

l_count$trait= factor(l_count$GENE, levels= unique(d$GENE))
l_count$supp= factor(l_count$supp, levels= c('Locus-level','Coloc'))

print('done')
p3= ggplot(l_count, aes(PP, GENE, fill= supp)) +
theme_cowplot(font_size= 8) +
geom_col(alpha= 0.7) +
geom_hline(yintercept= 0) +
scale_fill_manual(values= c(colorBlindBlack8[4], colorBlindBlack8[2]), guide= F) +
theme(	axis.line= element_blank(),
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	panel.background = element_blank(),
	panel.border = element_rect(colour= 'black', fill= NA, size=1),
	axis.text.y= element_blank(),
	axis.ticks.y= element_blank(),
	axis.title= element_blank(),
	plot.margin = unit(c(0, 0.1, 0, 0), "cm")) +
scale_x_continuous(limits= c(0, 10), expand= c(0,0), labels= seq(0,10, 2), breaks= seq(0, 10, 2))

x1= plot_grid(p1, p3, nrow= 1, align= 'h', rel_widths= c(2, 0.5))
x2= plot_grid(p1, p2, nrow= 2, align= 'v', rel_heights= c(2, 0.3))

ggsave(snakemake@output[[1]], plot= x1, width= 127 - 1, height= 127 - 25 - 1, units= 'mm', dpi= 300)
ggsave(snakemake@output[[2]], plot= x2, width= 103 - 1, height= 127 - 25 - 1, units= 'mm', dpi= 300)

################## Genetic correlations

d= fread(snakemake@input[[2]])

d= filter(d, grepl('GAraw', p1), !grepl('BW', p2), !grepl('male', p2))
#d$p1= 'Gestational duration (maternal)'
d$p1= 'Maternal'
x= fread(snakemake@input[[2]])

x= filter(x, grepl('GA_fetal', p1), !grepl('BW', p2), !grepl('male', p2))
#x$p1= 'Gestational duration (fetal)'
x$p1= 'Fetal'
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
                ifelse(trait== 'BW_fetal_effect', 'Fetal only',
                ifelse(trait== 'Preeclampsia', 'Pre-eclampsia',
                ifelse(trait== 'BW_maternal_effect', 'Maternal only',
                ifelse(trait== 'PCOS', 'Polycystic ovary syndrome', 'Endometriosis'))))))))))))))))))))))


d= filter(d, trait!= 'GA fetal effect')

d$cluster= with(d, ifelse(trait %in% pregnancy, 'Pregnancy', ifelse(trait %in% uterus, 'Reproductive organs', ifelse(trait %in% fitness, 'Fitness', 'Sex-hormone related'))))

d$colour= with(d, ifelse(cluster== 'Pregnancy', colorBlindBlack8[3], ifelse(cluster== 'Reproductive organs', colorBlindBlack8[1], ifelse(cluster== 'Fitness', colorBlindBlack8[7], colorBlindBlack8[8]))))

d= arrange(d, cluster)

d$trait= factor(d$trait, levels= traits)

colors <- filter(d, !duplicated(trait)) %>% arrange(trait) %>% pull(colour)


d$sig= ifelse(d$p< 0.05/ (nrow(d)/2), '**', ifelse(d$p< 0.05, '*', ''))
d= filter(d, p1== 'Maternal')
d$p1= 'Gestational duration'

rg_plot= ggplot(d, aes(trait, p1, fill= rg)) +
geom_tile(colour = "white", size= 1) +
theme_cowplot(font_size= 9) +
scale_fill_gradient2(low= colorBlindBlack8[2], high= colorBlindBlack8[4], mid= 'white', guide= F) +
theme(axis.text.x = element_text(angle = 45, hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
scale_x_discrete(position = "top") +
geom_text(data= d, aes(trait, p1, label= sig), size= 6/ .pt) +
theme(  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
	axis.ticks= element_blank(),
        panel.border = element_rect(colour= 'black', fill= NA, size=1),
        plot.margin = unit(c(0, 1, 0, 0), "cm"),
        axis.line= element_blank(),
	axis.text.x= element_text(angle= 45, hjust=0, colour= colors))


x2= plot_grid(rg_plot,p1, nrow= 2, align= 'v', rel_heights= c(0.85, 2))

ggsave(snakemake@output[[3]], plot= x2, width= 113 - 2.5, height= 127 - 25 - 1 , units= 'mm', dpi= 300)


