library(data.table)
library(dplyr)
library(ggplot2)
library(pROC)

df <- read.table(snakemake@input[[1]], header=T, sep="")
pheno <- read.table(snakemake@input[[2]], header=T, sep="\t")

df <- read.table('~/hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Data/Genotypes/best_beta_validation.profile', header=T, sep="")
pheno <- read.table('~/hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Training/Data/GAraw_pheno_validation.txt', header=T, sep="\t")

df <- subset(df, select=c("IID","SCORESUM"))
df <- merge(df, pheno, by="IID", all=T)

df <- df %>% rename("PRS"="SCORESUM")

df$preterm[df$GAraw<259] <-1
df$preterm[df$GAraw>273 & df$GAraw<=286] <-0

# create prs z score

mean(df$PRS)
# 0.1774635

se <- function(x) sd(x) / sqrt(length(x)) 
se(df$PRS)
#0.01671047

df$prs_Z <- (df$PRS - 0.1774635)/0.01671047

model1 <- lm(GAraw ~ prs_Z, data=df)

m1 <- summary(model1)$coefficients[1:2,]
r1 <- as.data.frame(summary(model1)$r.squared)
ar1 <- as.data.frame(summary(model1)$adj.r.squared)
fm1 <- cbind(m1,r1,ar1)

fm1 <- tibble::rownames_to_column(fm1, "Predictors")

model2 <- lm(GAraw ~ prs_Z + PC1 + PC2 + PC3 + PC4 + PC5, data=df)

m2 <- summary(model2)$coefficients[1:7,]
r2 <- as.data.frame(summary(model2)$r.squared)
ar2 <- as.data.frame(summary(model2)$adj.r.squared)
fm2 <- cbind(m2,r2,ar2)

fm2 <- tibble::rownames_to_column(fm2, "Predictors")

model3 <- lm(GAraw ~ prs_Z + Batch + PC1 + PC2 + PC3 + PC4 + PC5, data=df)

m3 <- summary(model3)$coefficients[1:10,]
r3 <- as.data.frame(summary(model3)$r.squared)
ar3 <- as.data.frame(summary(model3)$adj.r.squared)
fm3 <- cbind(m3,r3,ar3)

fm3 <- tibble::rownames_to_column(fm3, "Predictors")

model4 <- lm(GAraw ~ prs_Z + Batch + nullipara + PC1 + PC2 + PC3 + PC4 + PC5, data=df)
m4 <- summary(model4)$coefficients[1:11,]
r4 <- as.data.frame(summary(model4)$r.squared)
ar4 <- as.data.frame(summary(model4)$adj.r.squared)
fm4 <- cbind(m4,r4,ar4)

fm4 <- tibble::rownames_to_column(fm4, "Predictors")

model5 <- lm(GAraw ~ prs_Z*nullipara + Batch + PC1 + PC2 + PC3 + PC4 + PC5, data=df)

m5 <- summary(model5)$coefficients[1:12,]
r5 <- as.data.frame(summary(model5)$r.squared)
ar5 <- as.data.frame(summary(model5)$adj.r.squared)
fm5 <- cbind(m5,r5,ar5)

fm5 <- tibble::rownames_to_column(fm5, "Predictors")

dfn <- subset(df, nullipara==1)

dfm <- subset(df, nullipara==0)

model6 <- lm(GAraw ~ prs_Z + Batch + nullipara + PC1 + PC2 + PC3 + PC4 + PC5, data=dfn)

m6 <- summary(model6)$coefficients[1:10,]
r6 <- as.data.frame(summary(model6)$r.squared)
ar6 <- as.data.frame(summary(model6)$adj.r.squared)
fm6 <- cbind(m6,r6,ar6)

fm6 <- tibble::rownames_to_column(fm6, "Predictors")

model7 <- lm(GAraw ~ prs_Z + Batch + nullipara + PC1 + PC2 + PC3 + PC4 + PC5, data=dfm)

m7 <- summary(model7)$coefficients[1:10,]
r7 <- as.data.frame(summary(model7)$r.squared)
ar7 <- as.data.frame(summary(model7)$adj.r.squared)
fm7 <- cbind(m7,r7,ar7)

fm7 <- tibble::rownames_to_column(fm7, "Predictors")

##### Preterm birth
###################


model8 <- glm(preterm ~ prs_Z, data=df, family=binomial())

#### Get odds ratio and confidence intervals for model 8
m8 <- exp(cbind(OR = coef(model8), confint(model8)))



# Get AUC from model 8 
predicted <- predict(model8, df, type="response")
rocobj1 <- roc(df$preterm, predicted)

AUC <- (ci(rocobj1))

auc1 <- as.data.frame(rbind(m8, AUC))

m8 <- as.data.frame(summary(model1_1)$coefficients[1:2,])
me8 <- exp(cbind(coef(model8), confint(model8)))
m8 <- cbind(m8, me8)
m8 <- m8 %>% rename(OR = V1)

m8 <- tibble::rownames_to_column(m1_1, "Predictors")



write.table(auc1, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Results/preterm_model1.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

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
                    theme(text = element_text(size=20)) +
                    theme(
                      legend.position = c(.95, .5),
                      legend.justification = c("right", "top"),
                      legend.box.just = "right",
                      legend.background = element_rect(fill = "transparent"),
                      legend.text=element_text(size=18)) + 
                    ggplot2::annotate("text", 
                                    best.coords$specificity, 
                                    label = "0.57",
                                    color = "red",
                                    y= best.coords$specificity - 0.07,
                                    x = 0,
                                    size = 6) +
                    ggplot2::annotate("text", 
                                      best.coords$sensitivity, 
                                      label = "0.65",
                                      color = "blue",
                                      y= best.coords$sensitivity + 0.07,
                                      x = 0,
                                      size = 6) +
                    ggplot2::annotate("text", 
                                      best.coords$threshold, 
                                      label = "0.04",
                                      color = "black",
                                      x= best.coords$threshold + 0.025,
                                      y = -0.01,
                                      size=6)

png(file="./hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Plots/PTD_Sens_Spec.png",
    width=700, height=500)
senspec_plot
dev.off()


rm(AUC)

model9 <- glm(preterm ~ prs_Z + Batch + PC1 + PC2 + PC3 + PC4 + PC5, data=df, family=binomial())

#### Get odds ratio and confidence intervals for model 8
m9 <- exp(cbind(OR = coef(model9), confint(model9)))



#### Get AUC from model 8 need to complete this 
predicted2 <- predict(model9, df, type="response")
rocobj2 <- roc(df$preterm, predicted2)

AUC <- (ci(rocobj2))

auc2 <- as.data.frame(rbind(m9, AUC))

write.table(auc2, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Results/preterm_model2.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)



write.table(fm1, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Results/model1.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm2, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Results/model2.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm3, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Results/model3.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm4,'./hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Results/model4.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm5, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Results/model5.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm6, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Results/model6.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm7, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Results/model7.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(auc1, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Results/model8.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm9, './hunt-cloud/mnt/work2/chrisf/meta_gest_duration/prs_80-20_split/GAraw/Validation/Results/model9.txt',
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)





write.table(fm1, file = snakemake@output[[1]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm2, file = snakemake@output[[2]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm3, file = snakemake@output[[3]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm4, file = snakemake@output[[4]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm5, file = snakemake@output[[5]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm6, file = snakemake@output[[6]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm7, file = snakemake@output[[7]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(auc1, file = snakemake@output[[8]],
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(fm9, file = snakemake@output[[9]],
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

