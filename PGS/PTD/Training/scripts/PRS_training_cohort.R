

library(data.table)
library(dplyr)

#child <- fread('./hunt-cloud/mnt/archive/moba/pheno/v10/V10_1.1.1-200701/child.gz')
child <- fread(file = snakemake@input[[1]])

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

#kin <- fread('./hunt-cloud/mnt/archive/MOBAGENETICS/genotypes-base/aux/pedigree/mobagen-ethnic-core-samples.kin0')
kin= fread(file = snakemake@input[[2]])
kin = as.data.table(kin)

#flag <- fread('./hunt-cloud/mnt/archive/MOBAGENETICS/genotypes-base/aux/flaglist-merged/mobagen-flaglist-n99259.txt')
flag= fread(file = snakemake@input[[3]])
flag=as.data.table(flag)


#keep good genotypes and phenotypes

flag= filter(flag, genotypesOK== TRUE, phenoOK== TRUE)# I would also filter to keep only mothers here
flag <- merge(flag, relate, by="IID", all=FALSE)


####### Create Testing Dataframe 

###remove stillborns

child <- filter(child, (perinatal_stillborn=="No") & mother_genotyping_batch!="Ted") 

#del <- fread('./hunt-cloud/mnt/archive/moba/pheno/v10/V10_1.1.1-200701/delivery.gz')
del <- fread(file = snakemake@input[[4]],  sep="\t", header = TRUE)

del <- filter(del, (abruptio_placentae=="No" & placenta_previa=="No" & amniotic_polyhydramnion=="No" & amniotic_oligohydramnion=="No" &
                    (is.na(ceasarean)|ceasarean=="Emergency caesarean") & initiation!="Induced" & is.na(birth_defect) & is.na(art) &
                    position!="Breech" & mother_genotyping_batch!="Ted"))


#mh <- fread('./hunt-cloud/mnt/archive/moba/pheno/v10/V10_1.1.1-200701/mother_health.gz')
mh <- fread(file = snakemake@input[[5]],  sep="\t", header = TRUE)

mh <- filter(mh, (mother_chronic_hypertension=="No" & is.na(mother_diabetes) & mother_rheumatoid_arthritis=="No" & mother_genotyping_batch!="Ted"))


#preg <- fread('./hunt-cloud/mnt/archive/moba/pheno/v10/V10_1.1.1-200701/pregnancy.gz')
preg <- fread(file = snakemake@input[[6]],  sep="\t", header = TRUE)

preg <- filter(preg, (is.na(mother_smoking_end_cigarettes_per_day) & plural_birth=="Single birth" &
                        mother_pregnancy_hypertension=="No" & mother_pregnancy_early_preeclampsia=="No" & is.na(mother_pregnancy_preeclampsia) &
                        mother_pregnancy_eclampsia=="No" & mother_pregnancy_hellp=="No" & mother_genotyping_batch!="Ted"))


child1 <- child %>% select(child_SentrixID, mother_SentrixID, father_SentrixID, child_genotyping_batch, mother_genotyping_batch, father_genotyping_batch)
del1 <- del %>% select(child_SentrixID, mother_SentrixID, father_SentrixID, child_genotyping_batch, mother_genotyping_batch, father_genotyping_batch, birth_year)
mh1 <- mh %>% select(child_SentrixID, mother_SentrixID, father_SentrixID, child_genotyping_batch, mother_genotyping_batch, father_genotyping_batch)
preg1 <- preg %>% select(child_SentrixID, mother_SentrixID, father_SentrixID, child_genotyping_batch, mother_genotyping_batch, father_genotyping_batch, pregnancy_duration)

mum <- merge(child1, del1, all=FALSE)
mum <- merge(mum, mh1, all=FALSE)
mum <- merge(mum, preg1, all=FALSE)

mum <- mum[with(mum, order(birth_year)), ]

mum$dup <- duplicated(mum$mother_SentrixID)

mum <- filter(mum, (dup==FALSE))

mum <- mum %>% select(mother_SentrixID, mother_genotyping_batch, pregnancy_duration)
mum$IID <- mum$mother_SentrixID
mum <- mum %>% rename (
   "Batch" = "mother_genotyping_batch")

#good <- fread('./hunt-cloud/mnt/archive/MOBAGENETICS/genotypes-base/aux/flaglist-merged/mobagen-flaglist-n99259.txt')
good <- fread(file = snakemake@input[[3]],  sep="\t", header = TRUE)

goodmum <- merge(mum, good, by="IID", all=FALSE)
goodmum$FID <- goodmum$IID

###### create the phenotype files for PTD

goodmum <- goodmum %>% filter(!is.na(goodmum$pregnancy_duration))

goodmum$PTD[goodmum$pregnancy_duration<259] <-1
goodmum$PTD[goodmum$pregnancy_duration>273 & goodmum$pregnancy_duration<=286] <-0

goodmum <- goodmum %>% select(FID, IID, mother_SentrixID, Batch, PTD)

# remove missing (ie not preterm or full term)
goodmum <- na.omit(goodmum)

#pca <- fread('./hunt-cloud/mnt/archive/MOBAGENETICS/genotypes-base/aux/pca/mobagen-total/mobagen-total-proj-pc')
pca <- fread(file = snakemake@input[[7]],  sep="\t", header = TRUE)

pca <- pca %>% select(FID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

goodmum <- merge(goodmum, pca, by="FID", all=FALSE)


goodmum <- merge(goodmum, flag, by="IID", all=FALSE)


goodmum= goodmum %>% filter(!(IID %in% SelectRelated(kin, goodmum$IID)))


goodmum <- goodmum %>% select(FID, IID, mother_SentrixID, Batch, PTD, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)


#### Relatedness removed but duplicated mothers still in these IID - Remove

goodmum <- goodmum[!duplicated(goodmum$IID), ]


PTD <- subset(goodmum, select=c("FID","IID","PTD","Batch","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))

PTD$Batch= as.character(PTD$Batch)
PTD$Batch= ifelse(grepl('Norm', PTD$Batch), 'Norment', PTD$Batch)

print(table(PTD$Batch))


set.seed(101) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 80% of data as sample from total 'n' rows of the data
  
sample <- sample.int(n = nrow(PTD), size = floor(.8*nrow(PTD)), replace = F)
train <- PTD[sample, ]
validation  <- PTD[-sample, ]


print(nrow(train))
print(nrow(validation))


write.table(train, snakemake@output[[1]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(validation, snakemake@output[[2]], 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)















