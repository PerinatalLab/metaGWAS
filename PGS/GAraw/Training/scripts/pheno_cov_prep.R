

library(data.table)


phen <- read.table(file = snakemake@input[[1]],  sep="\t", header = TRUE)
rot2 <- read.table(file = snakemake@input[[2]],  sep="\t", header = TRUE)


phen <- merge(rot2, phen, by=c("FID", "IID"), all=FALSE)

phen$GAraw <- phen$pregnancy_duration

phen <- na.omit(phen)

phen$allPTD[phen$pregnancy_duration<259]<-1
phen$allPTD[phen$pregnancy_duration>258]<-0 

GAraw <- subset(phen, select=c("FID","IID","GAraw","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))

allPTD <- subset(phen, select=c("FID","IID","allPTD","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))

write.table(GAraw, snakemake@output[[1]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(allPTD, snakemake@output[[2]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
