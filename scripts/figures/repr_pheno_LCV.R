library(scales)
library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
library(tidyverse)
library(fmsb)

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

x= rbind(x, x1)

traits= unique(filter(x, p< 0.05/ 14, !grepl('BW', p2), !grepl('GA', p2)) %>% pull(p2))

d= fread(snakemake@input[[3]])

table_supp= d
table_supp$pheno= 'Gestational duration'
d$gcp.pm= ifelse(d$pval.gcpzero.2tailed< 0.05/length(traits), d$gcp.pm, 0)

d= filter(d, repr_pheno %in% traits)

d= arrange(d, desc(gcp.pm))

df= fread(snakemake@input[[4]])

table_supp2= df
table_supp2$pheno= 'Preterm delivery' 

table_supp= rbind(table_supp, table_supp2)

df$gcp.pm= ifelse(df$pval.gcpzero.2tailed< 0.05/length(traits), df$gcp.pm, 0)

df= filter(df, repr_pheno %in% traits)

d= inner_join(d, df, by= 'repr_pheno')
d$trait= d$repr_pheno
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

d$repr_pheno= d$trait
x= as.data.frame(matrix(d$gcp.pm.x, ncol= nrow(d)))
x=rbind(x, as.data.frame(matrix(d$gcp.pm.y, ncol= nrow(d))))



names(x)= d$repr_pheno
rownames(x)= c('Preterm delivery', 'Gestational duration ')
x= rbind(rep(1,nrow(d)) , rep(0,nrow(d)) , x)

inches= 25.4

pdf(snakemake@output[[1]], width= 88 / inches, height= 88 / inches)
par(mar=c(0,0,0,0))


radarchart(abs(x), axistype= 0,

    #custom polygon
    pcol= c(colorBlindBlack8[3], colorBlindBlack8[8]) , pfcol= c(alpha(colorBlindBlack8[3], 0.4), alpha(colorBlindBlack8[8], 0.4)) , plwd=1, pty= 16, plty= 1, vlcex= 0.8, vlabels= c('Testosterone\n(women)', 'Age at\nfirst birth', 'Age at\nmenopause', 'Number of\nlive births', 'SHBG\n(women)', 'CBAT\n(women)'),
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="#525252", caxislabels= seq(0, 1, 0.25), cglwd=0.8, calcex= 0.4

    #custom labels
    )
    
dev.off()

table_supp$trait= table_supp$repr_pheno
table_supp$trait= with(table_supp, ifelse(trait== 'GAraw', 'Maternal gestational duration',
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


fwrite(table_supp, snakemake@output[[2]], sep= '\t')
