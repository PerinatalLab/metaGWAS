library("dplyr")
library("knitr")
library("tidyr")
library(cowplot)
library(ggrepel)
library("data.table")
library('showtext')
library(ggtern)
options(warn=-1)


x= fread(snakemake@input[[1]], h= T)

colorBlindBlack8= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

x$rsid= with(x, ifelse(rsid== 'chrX:116013571', 'rs5991030', ifelse(rsid== 'chrX:132178061', 'rs5930554', rsid)))

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

d= fread(snakemake@input[[2]], header= T, select= c('RSID', 'ID', 'nearestGene'))

x= inner_join(x, d, by= c('rsid'= 'RSID'))

x$nearestGene= with(x, ifelse(rsid== 'rs3129768', 'HLA-DQA1', ifelse(rsid== 'rs5991030', 'AGTR2', ifelse(rsid== 'rs5930554', 'RAP2C', nearestGene)))) 

x$nearestGene= with(x, ifelse(rsid== 'rs6780427', 'KCNAB1', nearestGene))

x$nearestGene= with(x, ifelse(rsid== 'rs6879092', 'EBF1', nearestGene))

#d= gather(x, haplotype, beta, c('Paternal', 'MaternalT', 'MaternalNT'))

#d$rsid_label= with(d, paste0(RSID, ' (', nearestGene, ')'))


#max_beta= max(abs(d$beta))

#d$class= factor(d$class, levels= c("MF SD", "MF OD", "Maternal", "Fetal MatT", "Fetal"))

#d= arrange(d, class, desc(max_prob))

#d$rsid_label= factor(d$rsid_label, levels= unique(d$rsid_label))


colT= colorBlindBlack8[4]
colR= colorBlindBlack8[1]
colL= colorBlindBlack8[2]

x$MF= x$MF_OD + x$MF_SD
x$Fet= x$Fetal_MatT + x$Fetal

p1= ggtern(x, aes(Maternal, Fet, MF)) +
geom_point(colour= colorBlindBlack8[8], fill= colorBlindBlack8[8], shape= 21) +
scale_alpha_continuous(range= c(0.6, 1), guide= F) +
scale_size_continuous(range= c(.001, 10), guide= F) +
theme_custom(tern.plot.background = NULL, tern.panel.background = 'white', col.T = colT, col.L = colL, col.R = colR, col.grid.minor = "white") +
Tarrowlab("Fetal only effect") +
Larrowlab("Maternal only effect") +
Rarrowlab("Maternal and fetal effect")  +
theme_showarrows()  +
theme_notitles() +
theme(text=element_text(family="arial", size= 10),
        tern.axis.arrow.T = element_blank(),
        tern.axis.arrow.L = element_blank(),
        tern.axis.arrow.R = element_blank(),
        tern.axis.text.T = element_text(color = colT),
        tern.axis.text.L = element_text(color = colL),
        tern.axis.text.R = element_text(color = colR),
        tern.axis.arrow.text.T = element_text(color = colT),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        tern.axis.arrow.text.L = element_text(color = colL),
        tern.axis.arrow.text.R = element_text(color = colR),
        tern.panel.grid.major = element_line(linetype = 6, size = 0.3)) 
print('
ggtern(data=d, aes(-log10(pvalue_h1),-log10(pvalue_h2),-log10(pvalue_h3), label= nearestGene, size= abs(BETA), alpha= -log10(pvalue))) +
geom_point(colour= "black", fill= colorBlindBlack8[8], shape= 21) +
scale_alpha_continuous(range= c(0.6, 1), guide= F) +
scale_size_continuous(range= c(.001, 10), guide= F) +
theme_custom(tern.plot.background = NULL, tern.panel.background = "white", col.T = colT, col.L = colL, col.R = colR, col.grid.minor = "white") +
Tarrowlab("Maternal non-transmitted allele") + 
Larrowlab("Maternal transmitted allele") + 
Rarrowlab("Paternal transmitted allele")  +
theme_showarrows()  +
theme_notitles() +
theme(text=element_text(family="arial", size= 10),
	tern.axis.arrow.T = element_blank(),
	tern.axis.arrow.L = element_blank(),
	tern.axis.arrow.R = element_blank(),
        tern.axis.text.T = element_text(color = colT),
        tern.axis.text.L = element_text(color = colL),
        tern.axis.text.R = element_text(color = colR),
	tern.axis.arrow.text.T = element_text(color = colT), 
	plot.margin = margin(0, 0, 0, 0, "cm"), 
	tern.axis.arrow.text.L = element_text(color = colL),
	tern.axis.arrow.text.R = element_text(color = colR),
	tern.panel.grid.major = element_line(linetype = 6, size = 0.3)) +
geom_text(data= filter(d, nearestGene== "HAND2"), position= position_nudge_tern(y=0.05,x=-0.05/2,z=-0.05/2), aes(label=nearestGene), fontface= "bold", check_overlap=T, size= 8/ .pt, colour= "#525252", hjust= 1, vjust= 0.5)')


ggsave(snakemake@output[[1]], plot= p1, width= 95, height= 95, units= 'mm', dpi= 300)


d= select(x, rsid, ID, MF, Maternal, Fetal)

fwrite(x, snakemake@output[[2]], sep= '\t')
