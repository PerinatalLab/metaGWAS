
### Prepare workspace

library(data.table)
library(dplyr)


###########
##### Get the best betas and merge rsids
## Extract beta grid betas and merge with 


beta_grid <- read.table(file = snakemake@input[[1]],  sep="\t", header = TRUE)

r2_matrix <- read.table(file = snakemake@input[[2]],  sep="\t", header = TRUE)

result_grid <- read.table(file = snakemake@input[[3]],  sep="\t", header = TRUE)

info_snp <- read.table(file = snakemake@input[[4]],  sep="\t", header = TRUE)


##### Keep only the betas from the best model- compare the result_grid to the r2 matrix to determine which column is corresponding to the best model

match_grid <- which(r2_matrix$PRS_R2==result_grid$grid[1])


best_grid_beta <- subset(beta_grid, select=paste0('V',match_grid))

names(best_grid_beta)= "W_Beta"



##### Merge with the snps list (info_snps)

best_beta <- cbind(info_snp, best_grid_beta)

best_best <- best_beta %>% select(chr, pos, ID, rsid.ss, rsid, a0, a1, n_eff, EAF, beta, beta_se, p, W_Beta)

write.table(best_beta, snakemake@output[[1]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

df1 <- subset(best_beta, select=c("rsid", "a1"))

write.table(df1, snakemake@output[[2]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

df2 <- best_beta %>% select(rsid, W_Beta)

df2 <- inner_join(df1, df2)

write.table(df2, snakemake@output[[3]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)