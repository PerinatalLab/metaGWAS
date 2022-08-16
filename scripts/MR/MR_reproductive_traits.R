library(MendelianRandomization)
library(data.table)
library(dplyr)

if (!grepl('cluster', snakemake@output[[1]])){
d= fread(snakemake@input[[1]])
names(d)= c('ID', 'beta', 'se', 'pvalue', 'trait')
} else {
d= fread(snakemake@input[[1]])

}
x=fread(snakemake@input[[2]])
x= filter(x, !duplicated(ID))
d= inner_join(d, x, by= 'ID')



funk= function(temp_df){

inputMR= mr_input(bx = temp_df$beta,   bxse= temp_df$se,by = temp_df$BETA, byse = temp_df$SE)

if (nrow(temp_df)>3) {

z= mr_allmethods(inputMR)$Values
names(z)= c('method', 'estimate', 'se', 'lo95', 'up95', 'pvalue')
z$trait= unique(temp_df$trait)

} else {
z= mr_ivw(inputMR)

z= data.frame(method= 'IVW', estimate= z$Estimate, se= z$StdError, lo95= z$CILower, up95= z$CIUpper, pvalue= z$Pvalue, trait= unique(temp_df$trait))

}
return(z)
}


mr= lapply(split(d, d$trait), funk)

mr= do.call('rbind', mr)

fwrite(mr, snakemake@output[[1]], sep= '\t')

