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





d= fread('/mnt/hdd/common/pol/metaGWAS/repr_phenos/LDSC/results/repr_phenos_rg')


d$p2= gsub('.txt.sumstats.gz', '', apply(d[, 'p2'], 1, function(x) unlist(strsplit(x, 'LDSC/'))[2]))
d$p1= gsub('.txt.sumstats.gz', '', apply(d[, 'p1'], 1, function(x) unlist(strsplit(x, 'LDSC/'))[2]))

d$p1= d$p1
d$p1= with(d, ifelse(p1== 'miscarriage', 'Miscarriage',
                ifelse(p1== 'GA_fetal', 'GA fetal effect',
                ifelse(p1== 'BW_maternal', 'Birth weight maternal effect',
                ifelse(p1== 'AFB', 'Age at first birth',
                ifelse(p1== 'AMenarche', 'Age at menarche',
                ifelse(p1== 'AMenopause', 'Age at menopause',
                ifelse(p1== 'NLB', 'Number of live births',
                ifelse(p1== 'Testosterone_fem', 'Testosterone (women)',
                ifelse(p1== 'SHBG_fem', 'SHBG (women)',
                ifelse(p1== 'SHBG_male', 'SHBG (men)',
                ifelse(p1== 'CBAT_fem', 'CBAT (women)',
                ifelse(p1== 'CBAT_male', 'CBAT (men)',
                ifelse(p1== 'Oestradiol_fem', 'Oestradiol (women)',
                ifelse(p1== 'POP', 'Pelvic Organ Prolapse',
                ifelse(p1== 'Testosterone_male', 'Testosterone (men)',
                ifelse(p1== 'leiomyoma_uterus', 'Leiomyoma uterus',
                ifelse(p1== 'BW_fetal', 'Birth weight fetal effect',
                ifelse(p1== 'BW_fetal_effect', 'Birth weight fetal effect (adjusted MG)',
                ifelse(p1== 'Preeclampsia', 'Pre-eclampsia',
                ifelse(p1== 'BW_maternal_effect', 'Birth weight maternal effect (adjusted FG)',
                ifelse(p1== 'PCOS', 'Polycistic ovary syndrome', 'Endometriosis'))))))))))))))))))))))


d$trait= d$p2
d$trait= with(d, ifelse(trait== 'miscarriage', 'Miscarriage',
                ifelse(trait== 'GA_fetal', 'GA fetal effect',
                ifelse(trait== 'BW_maternal', 'Birth weight maternal effect',
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
                ifelse(trait== 'BW_fetal', 'Birth weight fetal effect',
                ifelse(trait== 'BW_fetal_effect', 'Birth weight fetal effect (adjusted MG)',
                ifelse(trait== 'Preeclampsia', 'Pre-eclampsia',
                ifelse(trait== 'BW_maternal_effect', 'Birth weight maternal effect (adjusted FG)',
                ifelse(trait== 'PCOS', 'Polycistic ovary syndrome', 'Endometriosis'))))))))))))))))))))))


ord= hclust( dist(d$rg, method = "euclidean"), method = "ward.D" )$order

d= d[ord, ]

d$p1= factor(d$p1, levels= unique(d$p1))

ggplot(d, aes(p1, p2, fill= rg)) + 
geom_tile() +
scale_fill_gradient2(low= colorBlindBlack8[2], high= colorBlindBlack8[4], mid= 'white', na.value = 'white') +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
	axis.title= element_blank(),
panel.background=element_rect(fill="white", colour="white"))

