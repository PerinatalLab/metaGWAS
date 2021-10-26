library(data.table)
library(parallel)
library(dplyr)
library(coloc)
library(tidyr)

top_pos= 123112292
dist_lim= 1.5*10**6

pph_outfile= snakemake@output[[1]]
results_outfile= snakemake@output[[2]]

cat('nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tpreg_trait\tpheno_FINNGEN\n', file = pph_outfile)
cat('snp\tV.df\tz.df1\tr.df1\tlABF.df1\tV.df2\tz.df2\tr.df2\tlABF.df2\tinternal.sum.lABF\tSNP.PP.H4\tpreg_trait\tpheno_FINNGEN\n', file= results_outfile)

prior1= 1 * 10**-4
prior2= 1 * 10**-4
prior12= 5 * 10**-6

format_FINNGEN= function(x){

temp_df= mani[mani$phenocode== x, ]

s_pheno= as.integer(temp_df[, 'n_cases']) / (as.integer(temp_df[, 'n_cases']) + as.integer(temp_df[, 'n_controls']))
pheno= as.character(unique(temp_df[, 'phenocode']))

N= (as.integer(temp_df[, 'n_cases']) + as.integer(temp_df[, 'n_controls']))
url= as.character(temp_df$path_https)

df= fread(url, select= c('#chrom', 'pos', 'ref', 'alt', 'beta', 'sebeta', 'maf', 'pval'))

names(df)= c('chr', 'pos', 'ref', 'alt', 'beta', 'se', 'maf', 'pvalue')

df= filter(df, pos %in% pos_hg38)

df$N= N
df$beta= ifelse(df$ref> df$alt, -1 * df$beta, df$beta)

df$ID_hg38= ifelse(df$ref> df$alt, paste(df$chr, df$pos, df$alt, df$ref, sep= ':'), paste(df$chr, df$pos, df$ref, df$alt, sep= ':'))

variants= data.frame(filter(df, pos== 123393445 | pos== 123346931))

if (nrow(variants)> 0){
variants$trait= ifelse(variants$pos== 123393445, 'Gestational duration', 'Birth weight, fetal effect')
}

variants$pheno= pheno


z= coloc_FINNGEN(d, df, 'Gestational duration', pheno, s_pheno)
z= coloc_FINNGEN(bw, df, 'Birth weight, fetal effect', pheno, s_pheno)

if (!file.exists(snakemake@output[[3]])){

fwrite(variants, snakemake@output[[3]], sep= '\t', row.names=F, col.names= T)

} else {
fwrite(variants, snakemake@output[[3]], sep= '\t', row.names=FALSE, col.names= FALSE, append= TRUE)

}


}

coloc_FINNGEN= function(trait1df, trait2df, trait, phenotype, prop_cases){

temp_df= inner_join(trait1df, trait2df, by= 'ID_hg38')
temp_df= filter(temp_df, SE> 0, se>0)

data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'quant', snp= temp_df$ID, MAF= temp_df$MAF)

data2= list(beta= temp_df$beta, varbeta= temp_df$se**2, N=temp_df$N, type= 'cc', snp= temp_df$ID, s= prop_cases, MAF= temp_df$maf)

print(phenotype)
myres= tryCatch({suppressWarnings(coloc.abf(data1, data2, p1= prior1, p2= prior2, p12= prior12))}, error= function(e) { return(0)}

)
        if (length(myres)==1 ) {
        PPH= data.frame(nsnps= 0, PP.H0.abf= 0, PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, preg_trait= trait, pheno_FINNGEN= phenotype)
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, , preg_trait= trait, pheno_FINNGEN= phenotype)
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')
        } else {
        PPH= data.frame(t(myres[[1]]))
        PPH$preg_trait= trait
	PPH$pheno_FINNGEN= phenotype
        if ((PPH$PP.H3.abf + PPH$PP.H4.abf) >= 0.1) {
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= myres[[2]]
        res$preg_trait= trait
        res$pheno_FINNGEN= phenotype
	if ((PPH$PP.H3.abf + PPH$PP.H4.abf) >= 0.1) {
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
	}
        } else {
        print('Not enough power')
        }
        }


}

d= fread(snakemake@input[[1]])

pos_hg38= as.numeric(separate(d, ID_hg38, into= c(NA, 'hg_pos38', NA, NA), sep= ':') %>% pull(hg_pos38))

bw= fread(snakemake@input[[2]])

mani= fread(snakemake@input[[3]])

format_FINNGEN(snakemake@wildcards[['FINNGEN_pheno']])

#z= do.call('rbind', df_list)

#fwrite(z, snakemake@output[[3]], sep= '\t', row.names=F, col.names= T)
