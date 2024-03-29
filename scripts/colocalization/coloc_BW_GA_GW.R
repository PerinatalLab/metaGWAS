library(data.table)
library(dplyr)
library(coloc)
library(parallel)

pph_outfile= snakemake@output[[1]]
results_outfile= snakemake@output[[2]]

cat('nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tlocus\n', file = snakemake@output[[1]])

cat('snp\tV.df\tz.df1\tr.df1\tlABF.df1\tV.df2\tz.df2\tr.df2\tlABF.df2\tinternal.sum.lABF\tSNP.PP.H4\tlocus\n', file= snakemake@output[[2]])


prior1= 1 * 10**-4
prior2= 1 * 10**-4
prior12= 5 * 10**-6

d= fread(snakemake@input[[1]], select= c('ID', 'CHR', 'POS', 'TOTALSAMPLESIZE', 'BETA', 'SE', 'pvalue', 'EAF'))
d$MAF= ifelse(d$EAF>0.5, 1 - d$EAF, d$EAF)

x= fread(snakemake@input[[2]], select= c('ID', 'N','BETA', 'SE', 'pvalue', 'EAF'))

x$MAF= ifelse(x$EAF>0.5, 1- x$EAF, x$EAF)

names(x)= c('ID', 'N', 'beta', 'se', 'p', 'eaf', 'maf')

d= inner_join(d, x, by= 'ID')

if (sum(is.na(d$eaf)) == nrow(d)) {
d$maf= d$MAF
} 

z= fread(snakemake@input[[3]])

z$CHR= as.numeric(gsub('chr', '', z$chr))

z$locus= 1:nrow(z)


funk= function(i) {
        row= z[i,]
	locus= paste0('locus_', i)
        temp_df= filter(d, CHR== as.integer(row[, 'CHR']), POS >= as.integer(row[, 'start']), POS<= as.integer(row[, 'stop']))
	
	if (nrow(temp_df)== 0) { 
	PPH= data.frame(nsnps= 0, PP.H0.abf= 0,PP.H1.abf= 0,  PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, locus= locus)
	pph_list[[i]]= PPH
	res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, locus= locus)
	res_list[[i]]= res
	fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, locus= locus)
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')

	} else {
	temp_df= filter(temp_df, SE>0, se>0)
	if (grepl('PCOS|miscarriage|POP|endometriosis|Preeclampsia|leiomyoma_uterus', snakemake@input[[2]])) {
        if (grepl('PCOS', snakemake@input[[2]])) {s_pheno=  (1184 + 670 + 157 +658 +984 + 485 + 462 )/ (1184 + 670 + 157 +658 +984 + 485 + 462 + 5799 + 1379 +2807 +6774 +2963+ 407 + 96172)}
        if (grepl('miscarriage', snakemake@input[[2]])) {s_pheno=   49996 / ( 174109 + 49996)}
        if (grepl('POP', snakemake@input[[2]])) {s_pheno= 7053 / (57407 + 7053) }
        if (grepl('endometriosis', snakemake@input[[2]])) {s_pheno= 1496 / (192678 + 1496 )}
        if (grepl('Preeclampsia', snakemake@input[[2]])){ s_pheno= 4630/ (4630 + 373345)}
        if (grepl('leiomyoma_uterus', snakemake@input[[2]])){ s_pheno= ( 14569) / (85792 + 14569)}
        if (grepl('allPTD', snakemake@input[[1]])) {
        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'cc', snp= temp_df$ID, MAF= temp_df$MAF, s= 0.067)
        } else if (grepl('postTerm', snakemake@input[[1]])) {
        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'cc', snp= temp_df$ID, MAF= temp_df$MAF, s= 0.122)
        } else {data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'quant', snp= temp_df$ID, MAF= temp_df$MAF) }

        data2= list(beta= temp_df$beta, varbeta= temp_df$se**2, N=temp_df$N, type= 'cc', snp= temp_df$ID, s= s_pheno, MAF= temp_df$maf)

        } else { 
	if (grepl('allPTD', snakemake@input[[1]])) {
        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'cc', snp= temp_df$ID, MAF= temp_df$MAF, s= 0.067)
        } else if (grepl('postTerm', snakemake@input[[1]])) {
        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'cc', snp= temp_df$ID, MAF= temp_df$MAF, s= 0.122)
        } else {data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'quant', snp= temp_df$ID, MAF= temp_df$MAF) }

        data2= list(beta= temp_df$beta, varbeta= temp_df$se**2, N=temp_df$N, type= 'quant', snp= temp_df$ID, MAF= temp_df$maf)
}
myres= tryCatch({suppressWarnings(coloc.abf(data1, data2, p1= prior1, p2= prior2, p12= prior12))}, error= function(e) { return(0)}
)
	if (length(myres)==1 ) { 
	PPH= data.frame(nsnps= 0, PP.H0.abf= 0, PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, locus= locus)
        pph_list[[i]]= PPH
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, locus= locus)
	res_list[[i]]= res
	fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')
	next 
	} else {
	PPH= data.frame(t(myres[[1]]))
        PPH$locus= locus
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= myres[[2]]
        res$locus= locus
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)

}
}
}

mclapply(1:nrow(z), funk, mc.cores= 3)
