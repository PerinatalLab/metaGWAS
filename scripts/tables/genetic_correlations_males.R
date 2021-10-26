library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

x= fread(snakemake@input[[1]])

x$p1= gsub('.txt.sumstats.gz', '', apply(x[, 'p1'], 1, function(x) unlist(strsplit(x, 'LDscore/'))[2]))
x$p2= gsub('.txt.sumstats.gz', '', apply(x[, 'p2'], 1, function(x) unlist(strsplit(x, 'LDSC/'))[2]))

x1= fread(snakemake@input[[2]])

x1$p1= gsub('.txt.sumstats.gz', '', apply(x1[, 'p1'], 1, function(x) unlist(strsplit(x, 'LDscore/'))[2]))
x1$p2= gsub('.txt.sumstats.gz', '', apply(x1[, 'p2'], 1, function(x) unlist(strsplit(x, 'LDSC/'))[2]))
d= rbind(x, x1)


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


d= filter(d, grepl('men', trait), !grepl('women', trait))





fwrite(d, snakemake@output[[1]], sep= '\t')

