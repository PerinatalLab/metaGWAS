## May have to run before starting R
#### export LC_ALL="en_US.UTF-8"

### Prepare workspace

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

#### Merge phenotype and covariate files together

library(data.table)
library(dplyr)

pheno <- fread(file = snakemake@input[[1]]) ### CHANGE TO fread
pheno= as.data.table(pheno)

##### Load Hapmap3 SNPS - restrict the analysis to HapMap3 only SNPS

#####info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))

info <- fread(file = snakemake@input[[2]]) # CHANGE TO fread
info = as.data.table(info)

#### Load and transform the summary stats

sumstats <- fread(file = snakemake@input[[3]]) # CHANGE TO fread
sumstats = as.data.table(sumstats)

# LDpred 2 require the header to follow the exact naming
sumstats <- sumstats %>% rename(
    chr = CHR,
    pos = POS,
    rsid = RSID,
    a1 = EFF,
    a0 = REF,
    beta_se = SE,
    p = pvalue,
    beta = BETA,
    n_eff = TOTALSAMPLESIZE
    )


#### remove "chr" from CHR column and change to integer

#sumstats$chr <- gsub("[a-zA-Z ]", "", sumstats$chr)
#sumstats$chr <- as.integer(sumstats$chr)

# Filter out hapmap SNPs

sumstats <- sumstats[sumstats$rsid %in% info$rsid,]

##### Calculate the LD Matrix

# Get maximum amount of cores
NCORES <- nb_cores()
# Open a temporary file
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
fam.order <- NULL
# preprocess the bed file (only need to do once for each data set)

print('Right before reading genotypes...')
#if(!file.exists(snakemake@params[[1]])){
      snp_readBed(snakemake@input[[4]])
      obj.bigSNP <- snp_attach(snakemake@params[[1]])
#         } else {obj.bigSNP <- snp_attach(snakemake@params[[1]])}

file.remove(snakemake@params[[1]])

#snp_readBed("file = snakemake@params[[1]]")


print('Right after reading genotypes..')
# now attach the genotype object
#obj.bigSNP <- snp_attach("file = snakemake@params[[2]].rds")
# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
# perform SNP matching
info_snp <- snp_match(sumstats, map)

write.table(info_snp, snakemake@output[[1]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes
# Rename the data structures
CHR <- map$chr
POS <- map$pos
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = "/mnt/work2/chrisf/omni_files/Bolt_omni")


#### split the ROT2 cohort into test and validation cohorts - Use Moba for real deal

#set.seed(1)
#ind.val <- sample(nrow(genotype), 1000)
#ind.test <- setdiff(rows_along(genotype), ind.val)
 

# calculate LD

for (chr in 1:23) {
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    print(paste('Chromosome', chr))
    # Calculate the LD
    corr0 <- snp_cor(
            genotype,
            ind.col = ind.chr2,
            ncores = NCORES,
            infos.pos = POS2[ind.chr2],
            size = 3 / 1000
        )
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}

print('Chromosomes read')
# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
        c("family.ID", "sample.ID"),
        c("FID", "IID"))

write.table(fam.order, snakemake@output[[2]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#### Perform LD score regression

print('Start LD_score regression')

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]


write.table(h2_est, snakemake@output[[3]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

print('Finished LD score regression')
#####Calculate the null R2

# Reformat the phenotype file such that y is of the same order as the 
# sample ordering in the genotype file
# y <- pheno[fam.order, on = c("FID", "IID")] (this didn't work but below does)
y <- pheno[order(match(pheno$IID, fam.order$IID)),]
# Calculate the null R2
# use glm for binary trait 
# (will also need the fmsb package to calculate the pseudo R2)
#library(fmsb)

null.model <- paste("PC", 1:10, sep = "", collapse = "+") %>%
    paste0("PTD ~ Batch +", .) %>%
    as.formula %>%
    lm(., data = y) %>%
    summary
null.r2 <- null.model$r.squared

print('Null model performed.')

####################
##### Grid model

# Prepare data for grid model
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
grid.param <- expand.grid(p = p_seq,
            h2 = h2_seq,
            sparse = c(FALSE, TRUE))

# Get adjusted beta from grid model
print('Getting adjusted beta...')

beta_grid <- snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)

print('Adjusted beta from grid model - done')
#### Obtain PRS


###### Using genome wide bed files

genotype <- obj.bigSNP$genotypes

print('Calculate PRS.')
# calculate PRS for all samples
ind.test <- 1:nrow(genotype)
pred_grid <- big_prodMat(genotype, 
                            beta_grid, 
                            ind.col = info_snp$`_NUM_ID_`)

print('Pred grid done')
#y1   <- pheno$PTD 

#Scores can only be done on 0 and 1 phenotypes

#grid.param$score <- big_univLogReg(as_FBM(pred_grid[ind.test, ]), y1[ind.test])$score

##### Get the final perfromance of the LDpred Models


##### Grid model

reg.formula <- paste("PC", 1:10, sep = "", collapse = "+") %>%
    paste0("PTD ~ PRS + Batch +", .) %>%
    as.formula
reg.dat <- y
max.r2 <- 0
r2_matrix <- grid.param
                   
print('Running grid model')
#for (i in 1:ncol(pred_grid)){
#    print(paste('Column', i))
#	reg.dat$PRS <- pred_grid[, i]
#    grid.model <- lm(reg.formula, data = reg.dat) %>%
#        summary  
#    r2_matrix[i,4] <- grid.model$r.squared 
#    r2_matrix[i,5] <- (grid.model$r.squared - null.r2)
##    r2_matrix[i,6] <- null.r2 
#    if(max.r2 < grid.model$r.squared){
#        max.r2 <- grid.model$r.squared 
#        model <- data.table(reg.dat$PRS)
#        }
#    }

for (i in 1:ncol(pred_grid)){
    print(paste('Column', i))
        reg.dat$PRS <- pred_grid[, i]
        if (sum(is.na(pred_grid[, i]), na.rm= T)== length(pred_grid[, i])){
        r2_matrix[i,4] <- NA
        r2_matrix[i,5] <- NA
    r2_matrix[i,6] <- null.r2
    } else {
    grid.model <- lm(reg.formula, data = reg.dat) %>%
        summary
    r2_matrix[i,4] <- grid.model$r.squared
    r2_matrix[i,5] <- (grid.model$r.squared - null.r2)
    r2_matrix[i,6] <- null.r2
    if(max.r2 < grid.model$r.squared){
        max.r2 <- grid.model$r.squared
        model <- data.table(reg.dat$PRS)
        }
    }
    }


print('Grid model done')

r2_matrix <- r2_matrix %>% rename (
    R2 = V4,
    PRS_R2 = V5,
    Null.R2 = V6)

result_grid <- data.table(
    grid = max.r2 - null.r2,
    null = null.r2
  )

print('Saving results...')

grid_prs <- model

write.table(beta_grid, snakemake@output[[4]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(grid.param, snakemake@output[[5]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(pred_grid, snakemake@output[[6]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(result_grid, snakemake@output[[7]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(grid_prs, snakemake@output[[8]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(r2_matrix, snakemake@output[[9]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(map, snakemake@output[[10]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#file.remove(snakemake@params[[1]])



