

library(data.table)
library(dplyr)
library(ggplot2)
#library(Rmisc)


df <- read.table('~/hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/Validation/Data/Genotypes/best_beta_validation.profile', header=T, sep="")
pheno <- read.table('~/hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/Training/Data/GAraw_pheno_validation.txt', header=T, sep="\t")

df <- subset(df, select=c("IID","SCORESUM"))
df <- merge(df, pheno, by="IID", all=T)

df1 <- df %>% rename("PRS"="SCORESUM")

### Identify the column number that corresponds to result grid (the optimium grid model)
  

grid <- cbind(fam_order, grid_prs)

grid <- subset(grid, select=c(FID, IID, V1))

df1 <- merge(df1, grid, by=c("FID","IID"), all=F)

df1 <- df1 %>% rename(prs=V1)


df1$prs_quant= ntile(df1$PRS, 10)


alpha <- 0.05


dfGA <- df1 %>% 
  group_by(prs_quant) %>% 
  summarize(mean = mean(GAraw),
            lower = mean(GAraw) - qt(1- alpha/2, (n() - 1))*sd(GAraw)/sqrt(n()),
            upper = mean(GAraw) + qt(1- alpha/2, (n() - 1))*sd(GAraw)/sqrt(n()))


gaplot <- ggplot(dfGA, aes(x=prs_quant, y=mean)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
  geom_point() + xlab("PGS Quantiles") +
  geom_smooth(method="lm",se=F, size=.8, formula = y ~ poly(x,3)) +
  ylab("Gestational Duration (Weeks)") +
  scale_x_discrete(limits=c(1:10))+ 
  scale_y_discrete(limits=c(277:285)) + 
  theme_bw() + theme(text = element_text(size=20/.pt))

gaplot <- gaplot + geom_smooth(method="lm",se=F, size=.8, formula = y ~ poly(x,3))

weeks <- c("39+4","39+5", "39+6","40+0",  "40+1", "40+2", "40+3", "40+4", "40+5")

gaplot <- gaplot + scale_y_continuous(labels=weeks, breaks=277:285, limits=c(277,285))
  

ggsave(file="~/hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/Validation/Plots/prs_quantiles.png",
    plot=gaplot, width=5.0, height=3.0, units="in", dpi=600)
 











