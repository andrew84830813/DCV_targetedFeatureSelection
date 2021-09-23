library(diffCompVarRcpp)
library(selEnergyPermR)
library(simplexDataAugmentation)
library(DiCoVarML)
library(caret)
library(dplyr)
library(compositions)
library(foreach)
library(parallel)
library(readr)
library(tidyr)
library(stringr)
library(glmnet) # glmnet



###-------------------------------------------------------*
metaData = data.frame(readr::read_tsv(file = "Data/NICUNEC.WGS.sample_details.tsv"))
otuData = data.frame(readr::read_tsv(file = "Data/NICUNEC.WGS.taxon_abundance.tsv"))
f_name = paste0("NEC_seed",sd,"_permute",permute_labels)


## Filter NEC samples for NEC  in the that occurs a least a week from now
metaData1 = metaData %>% 
  filter(Days.of.period.NEC.diagnosed<=-7 ) %>% 
  filter(Age.at.sample.collection..days.>0) 

## Control Samples
metaData2 = metaData %>% 
  filter(is.na(Days.of.period.NEC.diagnosed)) %>% 
  filter(Age.at.sample.collection..days. < max(metaData1$Age.at.sample.collection..days.)) 
# hist(metaData2$Age.at.sample.collection..days.)
# hist(metaData1$Age.at.sample.collection..days.)


## Merge Samples
metaData3 = metaData2
md1 = rbind(metaData1,metaData3)
md = data.frame(X1 = md1$X1,Status = md1$Necrotizing.enterocolitis)
table(md$Status)
table(md1$Days.of.period.NEC.diagnosed)
hist(md1$Age.at.sample.collection..days.)
samples =md1 %>% 
  group_by(Subject.ID,Necrotizing.enterocolitis) %>% 
  summarise(n = n())

#Pre Process OTU Data
otuData = subset(otuData,select = c("X1",metaData$X1))
## Sep
otuData1 = tidyr::separate(otuData,col = 1,
                           into = c("Kingdom","Phylum","Class","Order",
                                    "Family","Genus","Species","fhfhf"),
                           remove = F,sep = ";")

## retain taxa  with Genus level
otuData1 = otuData1[!is.na(otuData1$Genus),]
otuData1 = otuData1 %>% 
  select(Genus,starts_with("SRR"))
otuData1[is.na(otuData1)] = 0
otuData1 = tidyr::gather(otuData1,"SampleID","Counts",2:ncol(otuData1))
otuData1 = otuData1 %>% 
  dplyr::group_by(SampleID,Genus) %>% 
  dplyr::summarise(counts = sum(Counts))
otuData1 = tidyr::spread(otuData1,key = "Genus","counts",fill = 0)
colnames(otuData1)[1] = c("X1")


## Append Metadata
otuData1 = left_join(md,otuData1)
df = data.frame(Status = otuData1$Status,otuData1[,-2:-1],row.names = otuData1$X1)

## Write Data Output
exp = list(mbiome = df,metadata = md1)
save(exp ,file = "Output/microbiomeDB_dataAtLeast7DaysBeforeNEC.Rda")
write_csv(x = df,file = "Output/microbiomeDB_mbiomeAtLeast7DaysBeforeNEC.csv")
write_csv(x = md1,file = "Output/microbiomeDB_metadataAtLeast7DaysBeforeNEC.csv")






## Define Common folds
fold_matrix = matrix(rep(0,nrow(df)*20),nrow = nrow(df))
k_fold = 5
for(s in 1:20){
  set.seed(s)
  overll_folds = caret::createFolds(samples$Necrotizing.enterocolitis,k =k_fold,list = F)
  samples1 = samples
  samples1$fold = overll_folds
  samples1 = left_join(md1,samples1)
  fold_matrix[,s] = samples1$fold
}
fold_matrix = data.frame(fold_matrix)
write_csv(x = fold_matrix,file = "Output/commonFolds_NEC_crossValidation.csv")


