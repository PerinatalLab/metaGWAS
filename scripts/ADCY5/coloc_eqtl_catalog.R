library(data.table)
library(dplyr)
library(coloc)
library(parallel)

pph_outfile= snakemake@output[[1]]
results_outfile= snakemake@output[[2]]

cat('nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tpreg_trait\tgene\teqtl_data\n', file = snakemake@output[[1]])

cat('snp\tV.df\tz.df1\tr.df1\tlABF.df1\tV.df2\tz.df2\tr.df2\tlABF.df2\tinternal.sum.lABF\tSNP.PP.H4\tpreg_trait\tgene\teqtl_data\n', file= snakemake@output[[2]])

prior1= 1 * 10**-4
prior2= 1 * 10**-4
prior12= 5 * 10**-6

eqtl_coloc= function(temp_df, trait, gene, eqtl_data){
if (nrow(temp_df)== 0) {
        
        PPH= data.frame(nsnps= 0, PP.H0.abf= 0,PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, preg_trait= trait, geneid= gene, eqtl_data_n= eqtl_data)
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
	res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, preg_trait= trait, geneid= gene, eqtl_data_n= eqtl_data)
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)    
        print('next')
        } else {
        data1= list(beta= temp_df$BETA, varbeta= temp_df$SE**2, N=temp_df$TOTALSAMPLESIZE, type= 'quant', snp= temp_df$rsid, MAF= temp_df$MAF)
        data2= list(beta= temp_df$beta, varbeta= temp_df$se**2, N=temp_df$N, type= 'quant', snp= temp_df$rsid, MAF= temp_df$maf)
        myres= tryCatch({suppressWarnings(coloc.abf(data1, data2, p1= prior1, p2= prior2, p12= prior12))}, error= function(e) { return(0)}
)       
        if (length(myres)==1 ) {  
        PPH= data.frame(nsnps= 0, PP.H0.abf= 0, PP.H1.abf= 0, PP.H2.abf= 0, PP.H3.abf= 0, PP.H4.abf= 0, preg_trait= trait, geneid= gene, eqtl_data_n= eqtl_data)
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= data.frame(snp= 'none', V.df1= 0, z.df1= 0, r.df1= 0, lABF.df1= 0, V.df2= 0, z.df2= 0, r.df2= 0, lABF.df2= 0, internal.sum.lABF= 0, SNP.PP.H4= 0, preg_trait= trait, geneid= gene, eqtl_data_n= eqtl_data)
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        print('next')
        } else {
        PPH= data.frame(t(myres[[1]]))
	PPH$preg_trait= trait
        PPH$geneid= gene
	PPH$eqtl_data_n= eqtl_data
        if ((PPH$PP.H3.abf + PPH$PP.H4.abf) >= 0.01) {
        fwrite(PPH, pph_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        res= myres[[2]]
	res$preg_trait= trait
        res$geneid= gene
	res$eqtl_data_n= eqtl_data
        fwrite(res, results_outfile, sep= '\t', row.names=F, col.names= F, quote=F, append= T)
        } else {
        print('Not enough power')
        }

}
}
}

format_eqtl= function(temp_df){
	gene= unique(temp_df$gene_id)
	trait= unique(temp_df$preg_trait)
	eqtl_data= unique(temp_df$eqtl_data)
#	variants= filter(temp_df, pos== 123112292 | pos== 123065778)
	variants= filter(temp_df, position== 123393445 | position== 123346931)
#	variants$preg_trait= ifelse(variants$pos== 123112292, 'Gestational duration', 'Birth weight, fetal effect')
	variants$preg_trait= ifelse(variants$pos== 123393445, 'Gestational duration', 'Birth weight, fetal effect')
	variants$gene= gene
	variants$eqtl_data= eqtl_data
	temp_df = filter(temp_df, SE>0, se> 0)
	print(nrow(temp_df))
	eqtl_coloc(temp_df, trait, gene, eqtl_data)
	return(variants)
	
}


d= fread(snakemake@input[[1]], select= c('rsid', 'MAF', 'BETA', 'SE', 'TOTALSAMPLESIZE'))

bw= fread(snakemake@input[[2]], select= c('rsid', 'MAF', 'BETA', 'SE', 'TOTALSAMPLESIZE'))

df= fread(snakemake@input[[3]], h= T)

#names(df)= c('gene_id', 'rsid', 'maf', 'N', 'ref', 'alt', 'pvalue', 'beta', 'se') #c('molecular_trait_id', 'chromosome', 'position', 'ref', 'alt', 'variant', 'ma_samples', 'maf', 'pvalue', 'beta', 'se', 'type', 'ac', 'an', 'r2', 'molecular_trait_object_id', 'gene_id', 'median_tpm', 'rsid')

eqtl_data_name= gsub('.txt.gz', '', unlist(strsplit(snakemake@input[[3]], '/'))[11])
df$eqtl_data= eqtl_data_name

#d= inner_join(d, df, by= 'ID')
d= inner_join(d, df, by= 'rsid')
d$preg_trait= 'Gestational duration'

#bw= inner_join(bw, df, by= 'ID')
bw= inner_join(bw, df, by= 'rsid')
bw$preg_trait= 'Birth weight'


z= mclapply(split(d, d$gene_id), format_eqtl, mc.cores= 3)
z= do.call('rbind', z)
fwrite(z, snakemake@output[[3]], sep= '\t', row.names=F, col.names= T)

z= mclapply(split(bw, bw$gene_id), format_eqtl, mc.cores= 3)
z= do.call('rbind', z)
fwrite(z, snakemake@output[[3]], sep= '\t', row.names=F, col.names= F, append= T)
