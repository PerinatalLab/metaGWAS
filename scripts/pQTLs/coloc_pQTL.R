library(data.table)
library(dplyr)
library(coloc)
library(parallel)

df= fread(snakemake@input[[1]], select= c('CHR', 'POS', 'ID', 'BETA', 'SE', 'TOTALSAMPLESIZE', 'EAF'))

df= filter(df, !duplicated(ID))

df$MAF= ifelse(df$EAF>0.5, 1 - df$EAF, df$EAF)

x= fread(snakemake@input[[2]], select= c('CHR', 'POS', 'nearestGene'))
x= x[, c('CHR', 'POS', 'nearestGene')]
names(x)= c('CHR', 'pos2', 'nearestGene')

df= inner_join(df, x, by= 'CHR')

df= filter(df, POS>= pos2 - 1.5*10**6, POS< pos2 + 1.5*10**6)

z= fread(snakemake@input[[3]], select= c('chr', 'pos', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'TotalSampleSize'))

z$Allele1= toupper(z$Allele1)
z$Allele2= toupper(z$Allele2)

z$ID= with(z, ifelse(Allele1 > Allele2, paste(chr, pos, Allele2, Allele1, sep= ':'), paste(chr, pos, Allele1, Allele2, sep= ':')))

z$maf= ifelse(z$Freq1> 0.5, 1 - z$Freq1, z$Freq1)

z= select(z, ID, maf, Effect, StdErr, TotalSampleSize)

df= inner_join(df, z, by= 'ID')

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
	protein= unique(temp_df$nearestGene)
        if (nrow(temp_df)== 0) {

        PPH= data.frame(nsnps= 0, PP.H0.abf= 0,PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, protein= protein)
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, protein= protein)
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')

        } else {
	temp_df = filter(temp_df, SE>0, StdErr> 0)

	if (grepl('allPTD', snakemake@input[[1]])) {
        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'cc', snp= temp_df$ID, s= 0.067, MAF= temp_df$MAF)
        } else if (grepl('postTerm', snakemake@input[[1]])) {
        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'cc', snp= temp_df$ID, s= 0.122, MAF= temp_df$MAF)
        } else {data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N= temp_df$TOTALSAMPLESIZE, type= 'quant', snp= temp_df$ID, MAF= temp_df$MAF) }

        data2= list(beta= temp_df$Effect, varbeta= temp_df$StdErr**2, N=temp_df$TotalSampleSize, type= 'quant', snp= temp_df$ID, MAF= temp_df$maf)
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



mclapply(split(df, df$nearestGene), colocalization_eqtl, mc.cores= 3)

