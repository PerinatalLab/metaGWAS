library(dplyr)
library(data.table)

d= fread(snakemake@input[[1]])
d= filter(d, !is.na(Z))

x= fread(snakemake@input[[2]])
x= filter(x, !is.na(Z))

ld= fread(snakemake@input[[3]])

d= inner_join(d, x, by= 'SNP')
d= inner_join(d, ld, by= 'SNP')


source(snakemake@params[[1]])
setwd(snakemake@params[[2]])

LCV= RunLCV(d$L2, d$Z.y, d$Z.x, ldsc.intercept= 1, n.1= (d$N.y), n.2= (d$N.x))

cat('zscore\tpval.gcpzero.2tailed\tgcp.pm\tgcp.pse\trho.est\trho.err\tpval.fullycausal1\tpval.fullycausal2\th2.zscore1\th2.zscore2\tpheno\trepr_pheno\n', file = snakemake@output[[1]])

z= data.frame(zscore= LCV$zscore, pval.gcpzero.2tailed= LCV$pval.gcpzero.2tailed, gcp.pm= LCV$gcp.pm, gcp.pse= LCV$gcp.pse, rho.est= LCV$rho.est, rho.err= LCV$rho.err, pval.fullycausal1= LCV$pval.fullycausal[1],pval.fullycausal2= LCV$pval.fullycausal[2], h2.zscore1= LCV$h2.zscore[1], h2.zscore2= LCV$h2.zscore[2], pheno= snakemake@wildcards[['pheno']], repr_pheno= snakemake@wildcards[['repr_pheno']])

fwrite(z, snakemake@output[[1]], sep= '\t')
