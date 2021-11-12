
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




d= fread(snakemake@input[[1]])
d1= fread(snakemake@input[[2]])

d$beta= as.numeric(d$beta)
d$se= as.numeric(d$se)
d$pvalue= as.numeric(d$pvalue)

d1= filter(d1, PP.H4.abf> 0.75)
d= filter(d, pheno %in% d1$pheno_FINNGEN)

mani= fread(snakemake@input[[3]], select= c('phenocode', 'name'))
names(mani)= c('pheno', 'description')

d= inner_join(d, mani, by= 'pheno')

x= fread(snakemake@input[[4]])

x1= fread(snakemake@input[[5]])
x1= filter(x1, PP.H4.abf> 0.75)

x= filter(x, pheno %in% x1$pheno_PAN_UKBB)

mani=fread(snakemake@input[[6]], select= c('phenocode', 'trait_type', 'description'))
mani$pheno= paste(mani$trait_type, mani$phenocode, sep= '_')

x= inner_join(x, mani, by= 'pheno')

d$zscore= d$beta / d$se
x$zscore= x$beta / x$se

d= select(d, pheno, description, zscore, pvalue, trait)
x= select(x, pheno, description, zscore, pvalue, trait)
d= bind_rows(d, x)




d$zscore= ifelse(d$zscore> 10, 10, ifelse(d$zscore< -10, -10, d$zscore))

d$trait= ifelse(d$trait== 'Gestational duration', 'rs28654158 (gestational duration)', 'rs11708067 (birth weight)')

d$trait= factor(d$trait, levels= rev(c('rs28654158 (gestational duration)', 'rs11708067 (birth weight)')))

d$description= with(d, ifelse(grepl('Other diabetes', description), 'Other diabetes', description))

d$description= with(d, ifelse(grepl('Non-insulin-dep', description), 'Non-insulin dependent diabetes', description))
d$description= with(d, ifelse(grepl('Diabetes, varying def', description), 'Diabetes, wide', description))
d$description= with(d, ifelse(grepl('Intestinal adhesions', description), 'Intestinal adhesions', description))

d$description= with(d, ifelse(grepl('Type 2 diabetes, strict', description), 'Type 2 diabetes', description))

d$description= with(d, ifelse(grepl('Type 2 diabetes with other specified/multiple/unspecified complications', description), 'Type 2 diabetes with complications', description))

d$description= with(d, ifelse(grepl('Diabetes, insuline treatment', description), 'Diabetes, insuline treatment', description))

d$description= with(d, ifelse(grepl('Creatinine', description), 'Creatinine in urine', description))

ord <- hclust( dist(d$zscore, method = "euclidean"), method = "ward.D" )$order
d= d[ord, ]
d$description= factor(d$description, levels= unique(d$description))


p1= ggplot(d, aes(y= trait, x= description, fill= round(zscore), alpha= factor(as.numeric(pvalue< 5e-6)))) +
geom_tile(colour = "white", size= 1) +
theme_cowplot(font_size= 8) +
scale_alpha_discrete(guide=F, range= c(0.3, 1)) +
scale_fill_gradient2(low= colorBlindBlack8[3], high= colorBlindBlack8[8], mid= 'white', guide= F) +
theme(  axis.text.x= element_text(hjust= 1, angle= 45),
	axis.text.y= element_text(),
	axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin= unit(c(t= 0, r= 0, b= 0, l= 0), unit= 'cm'),
	axis.line.x = element_line(size = 0.3),
        axis.line.y = element_line(size = 0.3),
        axis.ticks= element_line(size= 0.3)) +
coord_equal() +
labs(x = NULL, y = NULL)




ggsave(snakemake@output[[1]], p1, height= 100, width= 127, units= 'mm', dpi= 300)

fwrite(d, snakemake@output[[2]], sep= '\t')


p1= ggplot(d, aes(description, trait, fill= round(zscore), alpha= factor(as.numeric(pvalue< 5e-6)))) +
geom_tile(colour = "white", size= 1) +
theme_cowplot(font_size= 8) +
scale_alpha_discrete(guide=F) +
scale_fill_gradient2(low= colorBlindBlack8[3], high= colorBlindBlack8[8], mid= 'white') +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin= margin(t= 0, r= 0, l= 0, b= 0, unit= 'pt')) +
scale_x_discrete(position = "top")


ggsave(snakemake@output[[3]], p1, height= 100, width= 140, units= 'mm', dpi= 300)


