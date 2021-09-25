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



#setwd("/nas/longleaf/home/andrew84/rProjects/DiCoVarFS_project")


# Load Helper Functions  ---------------------------------------------------------------
fnames = dir("Helper_Functions/")
for(f in fnames){
  source(paste0("Helper_Functions/",f))
}

# Setup Cluster ---------------------------------------------------------------

## Detect Cores
# clus <- parallel::makeCluster(10)
# doParallel::registerDoParallel(clus)




# Read External Args ---------------------------------------------------------------

args = c(2,0)
args = commandArgs(trailingOnly = TRUE)
sd = as.numeric(args[1]) # random seed selection
permute_labels = as.logical(as.numeric(args[2])) #should be 0(False) or 1(True)


## set random seed
set.seed(sd)




## load data
load("Output/LODO/crc_data.Rda")
df = obj$data
md = obj$md
combs = read_csv(file = "Output/LODO/halfSplits_LODO_crc.csv")
f_name = paste0("crcLODO_seed",sd,"_permute",permute_labels)





## View Class Labels
table(df$Status)
combs = data.frame(combs)

## get read of samples with reads <10e6
rs = rowSums(df[,-1])
bool = rs>10e6
df = df[bool,]
md = md[bool,]
table(md$dataset)


## Get folds
fld = if_else(md$dataset %in% combs[,sd],1,2)
k_fold = max(fld)

## Partition Data
allData = DiCoVarML::lodo_partition(df = data.frame(df),
                                    dataset_labels = fld,
                                    seed = sd)
benchmark = data.frame()



for(f in 1:k_fold){
  
  max_sparsity = .9
  
  ## Extract Test/Train Splilt
  ttData = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                            fold = f,permLabels = permute_labels,
                                            maxSparisty = max_sparsity,
                                            extractTelAbunance = F)
  ## Compute Total Parts
  number_parts = ncol(ttData$train_Data);number_parts
  nsamps = nrow(ttData$train_Data)
  table(ttData$y_test)
  classes = as.character(unique(ttData$y_train))
  
  
  # DCV --------------------------------------------------------------
  
  perc_totalParts2Keep = .75
  num_sets = 5
  
  base_dims = ncol(ttData$train_Data)
  max_parts = round(perc_totalParts2Keep*base_dims)
  sets = round(seq(1,max_parts,length.out = num_sets))[-1]
  
  
  # Run K-Fold Cross Validation ---------------------------------------------
  inner_perf = data.frame()
  
  
  ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")
  #ensemble = c("ranger","xgbTree","xgbLinear")
  max_sparsity = .9
  
  
  
  
  ## Train final Model
  ## Pre-Process
  trainx = data.frame(fastImputeZeroes(ttData$train_Data,impFactor = ttData$imp_factor))
  testx = data.frame(fastImputeZeroes(ttData$test_data,impFactor = ttData$imp_factor)) 
  
  ## compute log ratios
  lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = ttData$y_train,trainx))
  lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = ttData$y_test,testx))
  
  cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train, 
                                      includeInfoGain = T, nfolds = 1, numRepeats = 1, 
                                      rankOrder = F)
  
  ## Apply targted feature selection method
  tar_Features = 50
  
  tar_dcv = targeted_dcvSelection(trainx = trainx,
                                  testx = testx,
                                  dcv = cc.dcv$lrs,lrs.train = lrs.train,lrs.test = lrs.test,
                                  y_label = ttData$y_train,
                                  seed = sd,
                                  ensemble = ensemble,
                                  y_test = ttData$y_test,
                                  tarFeatures = tar_Features,
                                  ts.id = ttData$test_ids, 
                                  max_sparsity = max_sparsity
  )
  
  perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                    Dataset = f_name,Seed = sd,Fold = f,
                    Approach = tar_dcv$Performance$Approach,AUC = as.numeric(tar_dcv$Performance$AUC),
                    number_parts = tar_dcv$Performance$number_parts,number_ratios = tar_dcv$Performance$number_ratios ,
                    comp_time = tar_dcv$Performance$comp_time,
                    base_dims = ncol(trainx)
  )
  benchmark = rbind(benchmark,perf)
  
}


write_csv(x = benchmark,file = paste0("Results/",f_name,".csv"))

