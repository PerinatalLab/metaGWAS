library(data.table)
library(dplyr)
library(metafor)

funk= function(pheno) {

d_temp= d[d$outcome== pheno, ]

df_list= lapply(c('MT', 'MNT', 'PT'), function(i){


df_temp= d_temp[d_temp$haplotype== i, ]
print(nrow(d_temp))
res.FE= rma(yi= beta, sei= se,  data= df_temp, method= "FE")

df= data.frame(beta= res.FE$beta, se= res.FE$se, pvalue= res.FE$pval, lo95= res.FE$ci.lb, up95= res.FE$ci.ub, het_pvalue= res.FE$QEp, outcome= pheno, haplotype= i)

print(df)

return(df)

})

df= do.call('rbind', df_list)

return(df)

}


moba= fread(snakemake@input[[1]])
decode= fread(snakemake@input[[2]])
hunt= fread(snakemake@input[[3]])

d= rbind(moba, decode)
d= rbind(d, hunt)

df_list= lapply(unique(d$outcome), funk)

x= do.call('rbind', df_list)

df= group_by(d, haplotype, outcome) %>% summarize(n= sum(n))

x= inner_join(x, df, by= c('haplotype', 'outcome'))

fwrite(x, snakemake@output[[1]], sep= '\t')


