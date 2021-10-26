library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(tidyr)
library(showtext)
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

font_add("arial", "arial.ttf", bold= 'arial_bold.ttf')

showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)


d= fread(snakemake@input[[1]])
names(d)[8]= 'phenocode'
mani= fread(snakemake@input[[2]])

trait_list= c('biomarkers', 'continuous', 'icd10')
mani= mani[mani$trait_type %in% trait_list, ]

mani= filter(mani, saige_heritability_EUR> 0.01)
mani= mani[order(mani$saige_heritability_EUR, decreasing= TRUE), ]
mani= mani[!duplicated(mani$phenocode), ]

mani$phenocode= paste(mani$trait_type, mani$phenocode, sep= '_')
mani= mani[, c('phenocode', 'description')]
mani= mani[!duplicated(mani$description), ]

d= inner_join(d, mani[, c('description', 'phenocode')], by= 'phenocode')
d$cohort= 'UKBB'

x= fread(snakemake@input[[3]])
names(x)[8]= 'phenocode'
mani= fread(snakemake@input[[4]])
mani= mani[, c('phenocode', 'name')]
names(mani)= c('phenocode', 'description')
mani= mani[!duplicated(mani$description), ]

x= inner_join(x, mani, by= 'phenocode')
x$cohort= 'FINNGEN'

d= rbind(d, x)
d= d[order(d$PP.H4.abf, decreasing= F), ]
d= filter(d, PP.H4.abf> 0.01, PP.H4.abf + PP.H3.abf> 0.75)

d$preg_trait= factor(d$preg_trait)
empty_bar <- 6
to_add <- data.frame( matrix(NA, empty_bar*nlevels(d$preg_trait), ncol(d)) )
colnames(to_add) <- colnames(d)
to_add$preg_trait <- rep(levels(d$preg_trait), each=empty_bar)
d <- rbind(d, to_add)
d <- d %>% arrange(preg_trait)


d$id= seq(1, nrow(d))

label_data= d
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar
label_data$hjust<-ifelse( angle < -90, 1, 0)


label_data$angle<-ifelse(angle < -90, angle+180, angle)

#d$id= factor(d$id, levels= d$id[order(d$PP.H4.abf)])

base_data= d %>%
  group_by(preg_trait) %>%
  filter(is.na(PP.H4.abf)) %>%
  summarize(start=min(id), end=max(id) ) %>%
  rowwise() %>%
  mutate(title=mean(c(start, end)))

arc100= rep(1, 2)
arc75= rep(0.75, 2)
arc50= rep(0.50, 2)
arc25= rep(0.25, 2)

p1= ggplot(d, aes(as.factor(id), PP.H4.abf, fill= preg_trait, alpha= PP.H4.abf)) +
geom_bar(stat="identity", colour= NA) +
scale_alpha_continuous(range= c(0.4, 1), guide= F) +
geom_segment(data=base_data, aes(x = end, y = arc100, xend = start, yend = arc100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = end, y = arc75, xend = start, yend = arc75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = end, y = arc50, xend = start, yend = arc50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = end, y = arc25, xend = start, yend = arc25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  annotate("text", x = ((base_data$end[1] + base_data$start[1]) / 2), y = c((0.25 + 0.05) , (0.50 + 0.05), (0.75 + 0.05) , (1 + 0.05)), label = c("0.25", "0.50", "0.75", "1") , color="grey", size=3 , angle=0, fontface="bold", hjust= 0.5) +
   annotate("text", x = ((base_data$end[2] + base_data$start[2]) / 2), y = c((0.25 + 0.05) , (0.50 + 0.05), (0.75 + 0.05) , (1 + 0.05) ), label = c("0.25", "0.50", "0.75", "1") , color="grey", size=3, angle=15, fontface="bold", hjust=0.5) +
ylim(-0.2, 2) + # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
theme_cowplot() +
scale_fill_manual(values=colorBlindBlack8[c(2,4)], guide= F) +
scale_colour_manual(values=colorBlindBlack8[c(2,4)], guide= F) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")    ) +
  coord_polar(start = 0) +
geom_text(data=filter(label_data, PP.H4.abf> 0.75), aes(x= factor(id), y=PP.H4.abf + 0.01, label=description, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= filter(label_data, PP.H4.abf> 0.750)$angle, inherit.aes = FALSE) +
theme(panel.grid = element_blank(),
axis.title = element_blank(),
axis.text = element_blank(),
axis.ticks = element_blank())

p1= save_plot(snakemake@output[[1]], p1, base_width= 8, base_height= 8)

fwrite(d, snakemake@output[[2]], sep= '\t')
