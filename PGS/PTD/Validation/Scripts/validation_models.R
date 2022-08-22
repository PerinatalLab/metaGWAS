library(data.table)
library(dplyr)
library(ggplot2)
library(pROC)

df <- read.table(snakemake@input[[1]], header=T, sep="")
pheno <- read.table(snakemake@input[[2]], header=T, sep="\t")

df <- read.table('~/hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/PTD/Validation/Data/Genotypes/best_beta_validation.profile', header=T, sep="")
pheno <- read.table('~/hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/PTD/Validation/Data/PTD_pheno_Validation.txt', header=T, sep="\t")

df <- subset(df, select=c("IID","SCORESUM"))
df <- merge(df, pheno, by="IID", all=T)

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

m1 <- as.data.frame(summary(model1)$coefficients[1:2,])
m1e <- exp(cbind(coef(model1), confint(model1)))
m1 <- cbind(m1, m1e)
m1 <- m1 %>% rename(OR = V1)

m1 <- tibble::rownames_to_column(m1, "Predictors")

model1_1 <- glm(PTD1 ~ prs_Z, data=df, family="binomial")

m1_1 <- as.data.frame(summary(model1_1)$coefficients[1:2,])
m1e_1 <- exp(cbind(coef(model1_1), confint(model1_1)))
m1_1 <- cbind(m1_1, m1e_1)
m1_1 <- m1_1 %>% rename(OR = V1)

m1_1 <- tibble::rownames_to_column(m1_1, "Predictors")

model2 <- glm(PTD ~ prs_Z + PC1 + PC2 + PC3 + PC4 + PC5, data=df, family="binomial")

m2 <- as.data.frame(summary(model2)$coefficients[1:7,])
m2e <- exp(cbind(coef(model2), confint(model2)))
m2 <- cbind(m2, m2e)
m2 <- m2 %>% rename(OR = V1)

m2 <- tibble::rownames_to_column(m2, "Predictors")

model3 <- glm(PTD ~ prs_Z + Batch + PC1 + PC2 + PC3 + PC4 + PC5, data=df, family="binomial")

m3 <- as.data.frame(summary(model3)$coefficients[1:10,])
m3e <- exp(cbind(coef(model3), confint(model3)))
m3 <- cbind(m3, m3e)
m3 <- m3 %>% rename(OR = V1)

m3 <- tibble::rownames_to_column(m3, "Predictors")

##### Get ROC curves and AUC


#### Get AUC from model 1 

predicted <- predict(model1, df, type="response")
rocobj1 <- roc(df$PTD, predicted)

AUC <- as.data.frame(ci(rocobj1))

auc1 <- as.data.frame(t(AUC))
colnames(auc1) <- c('2.5%','AUC','97.5%')

# get sensitiviy and specificity and plot

mycoords <- coords(rocobj1, "all")
best.coords <- coords(rocobj1, "best", best.method="youden")

senspec_plot <- ggplot() + 
  geom_line(data=mycoords, aes(y = specificity, x=threshold, color = "Specificity")) + 
  geom_line(data=mycoords, aes(y = sensitivity, x=threshold, color = "Sensitivity")) +
  theme(legend.position="top",
        legend.title=element_blank()) +
  scale_color_manual(values=c('Sensitivity'='blue', 'Specificity'='red')) +
  geom_vline(xintercept=best.coords$threshold, lty=2, col="black") +
  geom_hline(yintercept=best.coords$specificity, lty=2, col="red") +
  geom_hline(yintercept=best.coords$sensitivity, lty=2, col="blue") +
  xlab("Cut-Off") + 
  ylab("Performance") +
  xlim(0, 0.2) +
  theme(text = element_text(size=20)) +
  theme(
    legend.position = c(.95, .5),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.background = element_rect(fill = "transparent"),
    legend.text=element_text(size=18)) + 
  ggplot2::annotate("text", 
                    best.coords$specificity, 
                    label = "0.55",
                    color = "red",
                    y= best.coords$specificity - 0.068,
                    x = 0.02,
                    size = 6) +
  ggplot2::annotate("text", 
                    best.coords$sensitivity, 
                    label = "0.61",
                    color = "blue",
                    y= best.coords$sensitivity + 0.068,
                    x = 0.02,
                    size = 6) +
  ggplot2::annotate("text", 
                    best.coords$threshold, 
                    label = "0.05",
                    color = "black",
                    x= best.coords$threshold + 0.01,
                    y = -0.015,
                    size=6)

png(file="./hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/PTD/Validation/Plots/PTD_Sens_Spec.png",
    width=700, height=500)
senspec_plot
dev.off()



#### Get AUC from model 2 

predicted <- predict(model2, df, type="response")
rocobj2 <- roc(df$PTD, predicted)

AUC <- as.data.frame(ci(rocobj2))

auc2 <- as.data.frame(t(AUC))
colnames(auc2) <- c('2.5%','AUC','97.5%')


#### Get AUC from model 3 

predicted <- predict(model3, df, type="response")
rocobj3 <- roc(df$PTD, predicted)

AUC <- as.data.frame(ci(rocobj3))

auc3 <- as.data.frame(t(AUC))
colnames(auc3) <- c('2.5%','AUC','97.5%')

aucall <- rbind(auc1,auc2,auc3)
aucall <- tibble::rownames_to_column(aucall, "model")
aucall <- aucall %>% mutate(model=recode(model, "ci(rocobj1)" = "PRS Only",
                                                "ci(rocobj2)" = "PRS + 5PCs",
                                                "ci(rocobj3)" = "PRS + PC + Batch"))



write.table(m1, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/PTD/Validation/Results/model1.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


write.table(m1_1, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/PTD/Validation/Results/model1_inverted.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(m2, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/PTD/Validation/Results/model2.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(m3, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/PTD/Validation/Results/model3.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(aucall, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/PTD/Validation/Results/AUC.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)




write.table(m1, file = snakemake@output[[1]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(m2, file = snakemake@output[[2]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(m3, file = snakemake@output[[3]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(aucall, file = snakemake@output[[4]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)



