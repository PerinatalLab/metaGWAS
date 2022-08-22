
library(data.table)
library(dplyr)
library(ggplot2)
library(pROC)
library(patchwork)


# Using PTD gwas results

# df <- read.table(snakemake@input[[1]], header=T, sep="")
# ptdpheno <- read.table(snakemake@input[[2]], header=T, sep="\t")

df <- read.table('~/hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/PTD/Validation/Data/Genotypes/best_beta_validation.profile', header=T, sep="")
ptdpheno <- read.table('~/hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/PTD/Validation/Data/PTD_pheno_Validation.txt', header=T, sep="\t")

df <- subset(df, select=c("IID","SCORESUM"))
df <- merge(df, ptdpheno, by="IID", all=T)

df <- df %>% rename("PRS"="SCORESUM")

mean(df$PRS)
# 0.504838

se <- function(x) sd(x) / sqrt(length(x)) 
se(df$PRS)
#0.002814106

df$prs_Z <- (df$PRS - 0.504838)/0.002814106

#for the inverted PTD
df$PTD1[df$PTD==0] <-1
df$PTD1[df$PTD==1] <-0

model1 <- glm(PTD ~ prs_Z, data=df, family="binomial")

#### Get AUC from model 1 

predicted <- predict(model1, df, type="response")
rocobj1 <- roc(df$PTD, predicted)

AUC <- as.data.frame(ci(rocobj1))

auc1 <- as.data.frame(t(AUC))
colnames(auc1) <- c('2.5%','AUC','97.5%')

# get sensitiviy and specificity and plot

ptd_mycoords <- coords(rocobj1, "all")
ptd_best.coords <- coords(rocobj1, "best", best.method="youden")

ptd_gwas_senspec_plot <- ggplot() + 
                            geom_line(data=ptd_mycoords, aes(y = specificity, x=threshold, color = "Specificity")) + 
                            geom_line(data=ptd_mycoords, aes(y = sensitivity, x=threshold, color = "Sensitivity")) +
                            theme(legend.position="top",
                                  legend.title=element_blank()) +
                            scale_color_manual(values=c('Sensitivity'='blue', 'Specificity'='red')) +
                            geom_vline(xintercept=ptd_best.coords$threshold, lty=2, col="black") +
                            geom_hline(yintercept=ptd_best.coords$specificity, lty=2, col="red") +
                            geom_hline(yintercept=ptd_best.coords$sensitivity, lty=2, col="blue") +
                            xlab("Cut-Off") + 
                            ylab("Performance") +
                            xlim(0, 0.11) +
                            theme_classic() +
                            theme(text = element_text(size=20)) +
                            theme(legend.position = c(.95, .5),
                                  legend.justification = c("right", "top"),
                                  legend.box.just = "right",
                                  legend.background = element_rect(fill = "transparent"),
                                  legend.text=element_text(size=18),
                                  legend.title=element_blank()) +  
                            ggplot2::annotate("text", 
                                              ptd_best.coords$specificity, 
                                              label = "0.55",
                                              color = "red",
                                              y= ptd_best.coords$specificity - 0.05,
                                              x = 0.032,
                                              size = 6) +
                            ggplot2::annotate("text", 
                                              ptd_best.coords$sensitivity, 
                                              label = "0.61",
                                              color = "blue",
                                              y= ptd_best.coords$sensitivity + 0.05,
                                              x = 0.032,
                                              size = 6) +
                            ggplot2::annotate("text", 
                                              ptd_best.coords$threshold, 
                                              label = "0.05",
                                              color = "black",
                                              x= 0.053,
                                              y = 0,
                                              size=6)



# Using GA gwas results


# df1 <- read.table(snakemake@input[[1]], header=T, sep="")
# pheno <- read.table(snakemake@input[[2]], header=T, sep="\t")

df1 <- read.table('~/hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Data/Genotypes/best_beta_validation.profile', header=T, sep="")
gapheno <- read.table('~/hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Training/Data/GAraw_pheno_validation.txt', header=T, sep="\t")

df1 <- subset(df1, select=c("IID","SCORESUM"))
df1 <- merge(df1, gapheno, by="IID", all=T)

df1 <- df1 %>% rename("PRS"="SCORESUM")

df1$preterm[df1$GAraw<259] <-1
df1$preterm[df1$GAraw>273 & df1$GAraw<=286] <-0

# create prs z score

mean(df1$PRS)
# 0.1774635

se <- function(x) sd(x) / sqrt(length(x)) 
se(df1$PRS)
#0.01671047

df1$prs_Z <- (df1$PRS - 0.1774635)/0.01671047

model8 <- glm(preterm ~ prs_Z, data=df1, family=binomial())

#### Get AUC from model 8 
predicted <- predict(model8, df1, type="response")
rocobj2 <- roc(df1$preterm, predicted)



ga_mycoords <- coords(rocobj2, "all")
ga_best.coords <- coords(rocobj2, "best", best.method="youden")

ga_senspec_plot <- ggplot() + 
                    geom_line(data=ga_mycoords, aes(y = specificity, x=threshold, color = "Specificity")) + 
                    geom_line(data=ga_mycoords, aes(y = sensitivity, x=threshold, color = "Sensitivity")) +
                    theme(legend.position="top",
                          legend.title=element_blank()) +
                    scale_color_manual(values=c('Sensitivity'='blue', 'Specificity'='red')) +
                    geom_vline(xintercept=ga_best.coords$threshold + 0.0005, lty=2, col="black") +
                    geom_hline(yintercept=ga_best.coords$specificity, lty=2, col="red") +
                    geom_hline(yintercept=ga_best.coords$sensitivity, lty=2, col="blue") +
                    xlab("Cut-Off") + 
                    ylab("Performance") +
                    xlim(0, 0.11) +
                    theme_classic() +
                    theme(text = element_text(size=20)) +
                    theme(legend.position = c(.95, .5),
                          legend.justification = c("right", "top"),
                          legend.box.just = "right",
                          legend.background = element_rect(fill = "transparent"),
                          legend.text=element_text(size=18),
                          legend.title=element_blank()) + 
                    ggplot2::annotate("text", 
                                      ga_best.coords$specificity, 
                                      label = "0.57",
                                      color = "red",
                                      y= ga_best.coords$specificity - 0.05,
                                      x = 0.025,
                                      size = 6) +
                    ggplot2::annotate("text", 
                                      ga_best.coords$sensitivity, 
                                      label = "0.65",
                                      color = "blue",
                                      y= ga_best.coords$sensitivity + 0.05,
                                      x = 0.025,
                                      size = 6) +
                    ggplot2::annotate("text", 
                                      ga_best.coords$threshold, 
                                      label = "0.04",
                                      color = "black",
                                      x= ga_best.coords$threshold + 0.009,
                                      y = 0,
                                      size=6)



ga_senspec_plot <- ga_senspec_plot + ggtitle("A") + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

ptd_gwas_senspec_plot <- ptd_gwas_senspec_plot + ggtitle("B") + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.75, "cm"))

sensspec <- ga_senspec_plot + ptd_gwas_senspec_plot + plot_layout(ncol=2,widths=c(2,2))

png(file="./hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/Combined_PTD_Plots/supp_fig10.png",
    width=1100, height=500)
sensspec
dev.off()

