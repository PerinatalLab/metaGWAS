library(data.table)
library(dplyr)
library(tidyr)

SelectRelated= function(kin_path, sample_list, df, var){
  kin= read.table(kin_path, h=T, comment.char = "", sep= '\t')
 kin= kin %>% filter(KINSHIP>0.125)
 kin= kin %>% filter(X.FID1 %in% sample_list & FID2 %in% sample_list)
 if (nrow(kin)== 0) {return()}
 kin= kin %>% mutate(ID1= paste(X.FID1,ID1, sep= ":"),
                      ID2= paste(FID2, ID2, sep= ":")) %>% select(ID1, ID2, KINSHIP)
  kin_temp= kin
  colnames(kin_temp)= c("ID2","ID1","KINSHIP")
  kin_temp= rbind(kin_temp, kin)
  kin_temp= kin_temp %>% add_count(ID1)
  kin_temp= kin_temp %>% add_count(ID2)
  kin_temp= arrange(kin_temp, n, nn)
  to_keep= list()

  for (i in 1:nrow(kin_temp)) {
    if (kin_temp[i,"ID1"] %in% unlist(kin_temp[0:i,"ID2"])) {
      kin_temp[i,"ID2"]= "X"
    }
    else
      to_keep[[i]] <- kin_temp[["ID1"]][i]
  }
  to_remove= kin_temp %>% filter(!(ID1 %in% unlist(to_keep))) %>% select(ID1)
  to_remove= to_remove[!duplicated(to_remove$ID1),]
  to_remove= to_remove %>% separate(ID1, c('FID','IID'), sep=":")

  return(unlist(to_remove[,1]))
}

Remove_dup= function(var){
df= pheno[order(pheno[,var], decreasing= T),]
df= df[!duplicated(df$IID),]
df= df %>% filter(!(IID %in% SelectRelated(paste0(kin_path, moms_kin), df$IID, df, 'IID')))
return(df)
}


pheno_vars= c('BARN_PID', 'MOR_PID','FAR_PID', 'SVLEN_UL_DG', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'FAAR', 'PARITY0', 'spont', 'VEKT', 'FAAR', 'GAMETOD', 'SVLEN_SM_DG', 'SVLEN_DG')


mfr= fread(snakemake@input[[1]])
ids= fread(snakemake@input[[4]])

mfr$PREG_ID= paste('PREG', (mfr$MOR_PID), (mfr$BARN_PID), (mfr$FAR_PID), sep= '_')

mfr= filter(mfr, PREG_ID %in% ids$PREG_ID)

pc= fread(snakemake@input[[3]])
names(pc)= c('BARN_PID', 'X2','X','X1','PC1','PC2','PC3','PC4','PC5','PC6','PC7', 'PC8','PC9','PC10')

mfr= mutate(mfr, PARITY0= as.numeric(PARITET_MFR==0),
                spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                is.na(INDUKSJON_PROSTAGLANDIN) & is.na(INDUKSJON_ANNET) &
                is.na(INDUKSJON_OXYTOCIN) & is.na(INDUKSJON_AMNIOTOMI)),
		GAMETOD= ifelse(is.na(SVLEN_UL_DG), 1, 0))

mfr$SVLEN_DG= ifelse(is.na(mfr$SVLEN_UL_DG), mfr$SVLEN_SM_DG, mfr$SVLEN_UL_DG)
mfr= filter(mfr, is.na(FLERFODSEL), DODKAT<6 | DODKAT>10, !is.na(SVLEN_DG))

pheno= inner_join(mfr, pc, by= 'BARN_PID')

pheno= pheno %>% filter(!(MOR_PID %in% SelectRelated(snakemake@input[[2]], pheno$MOR_PID, pheno, 'MOR_PID')))
pheno= pheno %>% filter(!(BARN_PID %in% SelectRelated(snakemake@input[[2]], pheno$BARN_PID, pheno, 'BARN_PID')))
pheno= pheno %>% filter(!(FAR_PID %in% SelectRelated(snakemake@input[[2]], pheno$FAR_PID, pheno, 'FAR_PID')))

pheno= filter(pheno, spont== 1)

pheno= select(pheno, all_of(pheno_vars))

fwrite(pheno, snakemake@output[[1]], sep= '\t')
