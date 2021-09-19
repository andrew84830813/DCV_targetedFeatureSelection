
.libPaths("C:/Users/andrew84/Documents/R/win-library/V4.1")

library(curatedMetagenomicData)
library(caret)
library(compositions)
library(dplyr)


conditions = combined_metadata %>% 
  dplyr::group_by(study_name,study_condition) %>% 
  dplyr::summarise(n = dplyr::n())


## Extract Data
## CRC
f_name = "YachidaS_2019"
f_name[2] = "FengQ_2015"
f_name[3] = "HanniganGD_2017"
f_name[4] = "YuJ_2015"
f_name[5] = "ZellerG_2014"
f_name[6] = "ThomasAM_2018a"
f_name[7] = "ThomasAM_2018b"
f_name[8] = "VogtmannE_2016"
f_name[9] = "ThomasAM_2019_c"
f_name[10] = "WirbelJ_2018"
f_name[11] = "GuptaA_2019"
for(f in f_name){
  expResults = getCuratedMetagenDataV3_0(fname = f)
  outFile = paste0("Output/CRC_exp/curatedMetaGenome_",f,".Rda")
  save(expResults,file = outFile)
}

## Benchmark
f_name[1] = "ZhuF_2020" ## binary
f_name[2] = "JieZ_2017" ## binary
f_name[3] = "VilaAV_2018" ##binary
f_name[4] = "RubelMA_2020" ## binary
f_name[5] = "HallAB_2017" ## binary IBD
f_name[6] = "NielsenHB_2014"
f_name[7] = "KarlssonFH_2013" #multiclass
for(f in f_name){
  expResults = getCuratedMetagenDataV3_0(fname = f)
  outFile = paste0("Output/curatedMetaGenome_",f,".Rda")
  save(expResults,file = outFile)
}



## standard 
f_name = "LiJ_2014"
expResults = getCuratedMetagenDataV3_0(fname = f_name)
outFile = paste0("Output/curatedMetaGenome_",f_name,".Rda")
save(expResults,file = outFile)



