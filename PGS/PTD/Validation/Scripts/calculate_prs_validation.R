
#### Validating the PRS 

### Prepare workspace

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

#### Merge phenotype and covariate files together

library(data.table)
library(dplyr)

 

betas <- read.table(file = snakemake@input[[1]],  sep="\t", header = TRUE)


beta1 <- subset(betas, select=c("rsid", "W_Beta"))
names(beta1) <- c("rsid", "beta")


beta2 <- as.data.frame(t(beta1))



betas.keep <- beta2
ind.keep <- subset(beta1, select=c("rsid"))

snp_readBed(snakemake@input[[2]])

# now attach the genotype object using the new set of genotypes from the validation cohort

obj.bigSNP <- snp_attach(snakemake@params[[1]])

# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")


# Assign the genotype to a variable for easier downstream analysis
G <- obj.bigSNP$genotypes


prs <-snp_PRS(
            G,
            betas.keep = beta1$beta,
                           )

fam.order <- as.data.table(obj.bigSNP$fam)

fam.order <- cbind(fam.order, prs)

prs_cohort <- subset(fam.order, select=c(family.ID, V1))
names(prs_cohort) = c("IID", "prs")


write.table(prs_cohort, file = snakemake@output[[1]],
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

file.remove(snakemake@params[[1]])
