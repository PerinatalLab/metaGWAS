library(data.table)
library(dplyr)
library(coloc)
library(parallel)

df= fread(snakemake@input[[1]], select= c('ID', 'CHR', 'POS', 'BETA', 'SE', 'TOTALSAMPLESIZE', 'EAF'))

df= filter(df, CHR== '3', POS>= 121612292, POS<= 124612292)
df$MAF= ifelse(df$EAF>0.5, 1 - df$EAF, df$EAF)

bw= fread(snakemake@input[[3]], select= c('ID', 'CHR', 'POS', 'BETA', 'SE', 'N', 'EAF'))
names(bw)= c('ID', 'CHR', 'POS', 'BETA', 'SE', 'TOTALSAMPLESIZE', 'EAF')
bw= filter(bw, CHR== 3, POS>= 121612292, POS<= 124612292)
bw$MAF= ifelse(bw$EAF>0.5, 1 - bw$EAF, bw$EAF)


z= fread(snakemake@input[[2]])
z$n= 716
df= inner_join(df, z, by= 'ID')

bw= inner_join(bw, z, by= 'ID')

df$trait= 'Gestational duration'
bw$trait= 'Birth weight'

rm(z)

pph_outfile= snakemake@output[[1]]
results_outfile= snakemake@output[[2]]


cat('nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tgene\n', file = snakemake@output[[1]])

cat('snp\tV.df\tz.df1\tr.df1\tlABF.df1\tV.df2\tz.df2\tr.df2\tlABF.df2\tinternal.sum.lABF\tSNP.PP.H4\tgene\n', file= snakemake@output[[2]])

prior1= 1 * 10**-4
prior2= 1 * 10**-4
prior12= 5 * 10**-6


df= data.frame(df)

colocalization_eqtl= function(temp_df){
	protein= unique(temp_df$gene)
	trait= unique(temp_df$trait)
        if (nrow(temp_df)== 0) {

        PPH= data.frame(nsnps= 0, PP.H0.abf= 0,PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, preg_trait= trait, protein= protein)
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, preg_trait= trait, protein= protein)
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')


        } else {
	temp_df = filter(temp_df, SE>0, se> 0)
        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'quant', snp= temp_df$ID)
        data2= list(beta= temp_df$beta, varbeta= temp_df$se**2, N=temp_df$n, type= 'quant', snp= temp_df$ID)
        myres= tryCatch({suppressWarnings(coloc.abf(data1, data2, MAF=temp_df$MAF, p1= prior1, p2= prior2, p12= prior12))}, error= function(e) { return(0)}
)
        if (length(myres)==1 ) { 
        PPH= data.frame(nsnps= 0, PP.H0.abf= 0, PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, preg_trait= trait, protein= protein)
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, preg_trait= trait, protein= protein)
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')
        } else {
        PPH= data.frame(t(myres[[1]]))
	PPH$trait= trait
        PPH$protein= protein
        if ((PPH$PP.H3.abf + PPH$PP.H4.abf) >= 0.8) {
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= myres[[2]]
	res$trait= trait
        res$protein= protein
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        } else {
        print('Not enough power')
        }
}
}
}



mclapply(split(df, df$gene), colocalization_eqtl, mc.cores= 3)

mclapply(split(bw, bw$gene), colocalization_eqtl, mc.cores= 3)

z= filter(df, POS== 123112292 | POS== 123065778)
z1= filter(bw, POS== 123112292 | POS== 123065778)

z$trait= 'Gestational duration'
z1$trait= 'Birth weight'

z= rbind(z, z1)

fwrite(z, snakemake@output[[3]], sep= '\t')

