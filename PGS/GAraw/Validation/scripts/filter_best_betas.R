
library(data.table)
library(dplyr)


betas <- read.table(file=snakemake@input[[1]],  sep="\t", header = TRUE)

testbim <- read.table(file=snakemake@input[[2]], header=F, sep="\t")

testbim <- testbim %>% rename ("rsid"="V2")

testbim <- subset(testbim, select=rsid)

rtest <- merge(testbim, betas, by="rsid", all=F)

write.table(rtest, file = snakemake@output[[1]],
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)







