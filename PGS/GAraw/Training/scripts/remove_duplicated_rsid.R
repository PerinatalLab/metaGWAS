

library(data.table)
library(dplyr)


df <- read.table(file = snakemake@input[[1]],  sep="\t", header = TRUE)

df <- subset(df, RSID!="")

df$dup <- duplicated(df$RSID)
df$dup1 <- duplicated(df$ID)

df1 <- subset(df, dup=="FALSE")
df1 <- subset(df1, dup1=="FALSE")

write.table(df1, snakemake@output[[1]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
