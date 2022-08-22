

library(data.table)
library(dplyr)


child <- read.table(file = snakemake@input[[1]],  sep="\t", header = TRUE)

relate <- subset(child, !is.na(mother_SentrixID), select=c("mother_SentrixID"))

relate <- relate %>% rename (IID=mother_SentrixID)

#### Create a file of IID of non-related mothers

SelectRelated= function(kin, sample_list){
 kin= kin %>% filter(Kinship>0.125) # You can increase this number to 0.125 as suggest by LDpred
 kin= kin %>% filter(ID1 %in% sample_list, ID2 %in% sample_list) # Keep only samples that are in the file


if (nrow(kin)>0){
 kin= kin %>% select(ID1, ID2, Kinship)
  kin_temp= kin
  colnames(kin_temp)= c("ID2","ID1","Kinship")
  kin_temp= rbind(kin_temp, kin)
  kin_temp= kin_temp %>% add_count(ID1)
  kin_temp= kin_temp %>% add_count(ID2)
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
  #to_remove= to_remove %>% separate(ID1, c('FID','ID'), sep=":")
  return(unlist(to_remove[,1]))
}
}


kin= read.table(file = snakemake@input[[2]],  sep="\t", header = TRUE)
kin = as.data.table(kin)

flag= read.table(file = snakemake@input[[3]],  sep="\t", header = TRUE)
flag=as.data.table(flag)

#flag= fread(flag)

#kin= fread(kin)


flag= filter(flag, genotypesOK== TRUE, phenoOK== TRUE)# I would also filter to keep only mothers here
flag <- merge(flag, relate, by="IID", all=FALSE)




####### Create Testing Dataframe 

###remove stillborns

child <- filter(child, (perinatal_stillborn=="No") & mother_genotyping_batch!="Ted") 




del <- read.table(file = snakemake@input[[4]],  sep="\t", header = TRUE)

del <- filter(del, (abruptio_placentae=="No" & placenta_previa=="No" & amniotic_polyhydramnion=="No" & amniotic_oligohydramnion=="No" &
                    (is.na(ceasarean)|ceasarean=="Emergency caesarean") & initiation!="Induced" & is.na(birth_defect) & is.na(art) &
                    position!="Breech" & mother_genotyping_batch!="Ted"))


mh <- read.table(file = snakemake@input[[5]],  sep="\t", header = TRUE)

mh <- filter(mh, (mother_chronic_hypertension=="No" & is.na(mother_diabetes) & mother_rheumatoid_arthritis=="No" & mother_genotyping_batch!="Ted"))



preg <- read.table(file = snakemake@input[[6]],  sep="\t", header = TRUE)

preg <- filter(preg, (is.na(mother_smoking_end_cigarettes_per_day) & plural_birth=="Single birth" &
                        mother_pregnancy_hypertension=="No" & mother_pregnancy_early_preeclampsia=="No" & is.na(mother_pregnancy_preeclampsia) &
                        mother_pregnancy_eclampsia=="No" & mother_pregnancy_hellp=="No" & mother_genotyping_batch!="Ted"))


child1 <- subset(child, select=c(1:12))
del1 <- subset(del, select=c(1:12,15))
mh1 <- subset(mh, select=c(1:12))
preg1 <- subset(preg, select=c(1:12,19))

mum <- merge(child1, del1, all=FALSE)
mum <- merge(mum, mh1, all=FALSE)
mum <- merge(mum, preg1, all=FALSE)

mum <- mum[with(mum, order(birth_year)), ]

mum$dup <- duplicated(mum$mother_SentrixID)

mum <- filter(mum, (dup==FALSE))

mum <- subset(mum, select=c("mother_SentrixID", "mother_genotyping_batch", "pregnancy_duration"))
mum$IID <- mum$mother_SentrixID
mum <- mum %>% rename (
   "Batch" = "mother_genotyping_batch")

good <- read.table(file = snakemake@input[[3]],  sep="\t", header = TRUE)

goodmum <- merge(mum, good, by="IID", all=FALSE)
goodmum$FID <- goodmum$IID

goodmum <- subset(goodmum, select=c(10,1,3,4))

pca <- read.table(file = snakemake@input[[7]],  sep="\t", header = TRUE)

pca <- subset(pca, select = c(FID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

goodmum <- merge(goodmum, pca, by="FID", all=FALSE)

goodmum <- subset(goodmum, !is.na(goodmum$pregnancy_duration))

goodmum <- merge(goodmum, flag, all=FALSE)


goodmum= goodmum %>% filter(!(IID %in% SelectRelated(kin, goodmum$IID)))
goodmum$FID= goodmum$IID



goodmum <- subset(goodmum, select=c(1:16))


#### Relatedness removed but duplicated mothers still in these IID - Remove

goodmum <- goodmum[!duplicated(goodmum$IID), ]

###### create the phenotype files for GAraw and allPTD


goodmum$GAraw <- goodmum$pregnancy_duration

goodmum <- na.omit(goodmum)


GAraw <- subset(goodmum, select=c("FID","IID","GAraw","Batch","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))

GAraw$Batch= as.character(GAraw$Batch)
GAraw$Batch= ifelse(grepl('Norm', GAraw$Batch), 'Norment', GAraw$Batch)

print(table(GAraw$Batch))

#GAraw= filter(GAraw, Batch== snakemake@wildcards[['batch']])

set.seed(101) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 80% of data as sample from total 'n' rows of the data
  
sample <- sample.int(n = nrow(GAraw), size = floor(.8*nrow(GAraw)), replace = F)
train <- GAraw[sample, ]
validation  <- GAraw[-sample, ]


print(nrow(train))
print(nrow(validation))


write.table(train, snakemake@output[[1]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(validation, snakemake@output[[2]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)















