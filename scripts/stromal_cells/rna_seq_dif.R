library(data.table)
library(dplyr)
library(DESeq2)
library(tidyverse)

df_list= list()

flist= list.files(snakemake@params[[1]], 'CL', full.names=T)

for (i in 1:length(flist)){
d= fread(flist[i])
cname= unlist(strsplit(flist[i], '/'))[10]
d= select(d, Name, NumReads)

names(d)= c('Name', cname)
df_list[[i]]= d

}

x= df_list %>% reduce(left_join, by = "Name")

cols= data.frame(row.names= colnames(x)[2:7], condition= colnames(x)[2:7], subject= colnames(x)[2:7])

cols$condition= gsub('.txt', '', sapply(strsplit(cols$condition, '-'), tail, 1))
cols$subject= sapply(strsplit(cols$subject, '-'), head, 1)
cts= as.matrix(x[, 2:7])
row.names(cts)= x$Name

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = cols,
                              design= ~ subject + condition)

dds= DESeq(dds)

res= results(dds, name="condition_unt_vs_dec")

res= data.frame(res)
res$geneid= row.names(res)

fwrite(res, snakemake@output[[1]], sep= '\t')
