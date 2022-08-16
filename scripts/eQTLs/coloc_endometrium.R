library(data.table)
library(dplyr)
library(coloc)
library(parallel)

df= fread(snakemake@input[[1]], select= c('RSID', 'BETA', 'SE', 'TOTALSAMPLESIZE', 'EAF'))

df= filter(df, !duplicated(RSID))

df$MAF= ifelse(df$EAF>0.5, 1 - df$EAF, df$EAF)

z= fread(snakemake@input[[2]])
z$n= 206
z$maf= ifelse(z$Freq< 0.5, 1 - z$Freq, z$Freq)
df= inner_join(df, z, by= c('RSID'= 'SNP'))

rm(z)

pph_outfile= snakemake@output[[1]]
results_outfile= snakemake@output[[2]]


cat('nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tprotein\n', file = snakemake@output[[1]])

cat('snp\tV.df\tz.df1\tr.df1\tlABF.df1\tV.df2\tz.df2\tr.df2\tlABF.df2\tinternal.sum.lABF\tSNP.PP.H4\tprotein\n', file= snakemake@output[[2]])

prior1= 1 * 10**-4
prior2= 1 * 10**-4
prior12= 5 * 10**-6


df= data.frame(df)

colocalization_eqtl= function(temp_df){
	protein= unique(temp_df$Gene)
        if (nrow(temp_df)== 0) {

        PPH= data.frame(nsnps= 0, PP.H0.abf= 0,PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, protein= protein)
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, protein= protein)
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')

        } else {
	temp_df = filter(temp_df, SE>0, se> 0)

	if (grepl('allPTD', snakemake@input[[1]])) {
        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'cc', snp= temp_df$RSID, s= 0.067, MAF= temp_df$MAF)
        } else if (grepl('postTerm', snakemake@input[[1]])) {
        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'cc', snp= temp_df$RSID, s= 0.122, MAF= temp_df$MAF)
        } else {data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N= temp_df$TOTALSAMPLESIZE, type= 'quant', snp= temp_df$RSID, MAF= temp_df$MAF) }

        data2= list(beta= temp_df$b, varbeta= temp_df$se**2, N=temp_df$n, type= 'quant', snp= temp_df$RSID, MAF= temp_df$maf)
        myres= tryCatch({suppressWarnings(coloc.abf(data1, data2, p1= prior1, p2= prior2, p12= prior12))}, error= function(e) { return(0)}
)
        if (length(myres)==1 ) { 
        PPH= data.frame(nsnps= 0, PP.H0.abf= 0, PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, protein= protein)
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, protein= protein)
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')
        } else {
        PPH= data.frame(t(myres[[1]]))
        PPH$protein= protein
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= myres[[2]]
        res$protein= protein
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
}
}
}



mclapply(split(df, df$Gene), colocalization_eqtl, mc.cores= 3)

