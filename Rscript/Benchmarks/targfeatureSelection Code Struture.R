
perc_totalParts2Keep = .75
num_sets = 5

base_dims = ncol(ttData$train_Data)
max_parts = round(perc_totalParts2Keep*base_dims)
sets = round(seq(1,max_parts,length.out = num_sets))[-1]


# Run K-Fold Cross Validation ---------------------------------------------
inner_perf = data.frame()


ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")
ensemble = c("ranger","xgbTree","xgbLinear")
max_sparsity = .9

for(sd1 in 1:1){
  set.seed(sd1)
  k_fold = 2
  overll_folds = caret::createFolds(ttData$y_train,k = k_fold,list = F)
  innerfold_data = lodo_partition(data.frame(Status = ttData$y_train,ttData$train_Data),
                                  dataset_labels = overll_folds,
                                  sd1)
  
  
  ## Get within fold cross validated performance 
  
  for(f in 1:k_fold){
    
    ## Partition inner fold
    innerFold = DiCoVarML::extractTrainTestSplit(foldDataList = innerfold_data,
                                                 fold = f,
                                                 maxSparisty = max_sparsity,
                                                 extractTelAbunance = F)
    
    suppressMessages(suppressWarnings({
    
      ## Pre-Process
      trainx = data.frame(fastImputeZeroes(innerFold$train_Data,impFactor = innerFold$imp_factor))
      testx = data.frame(fastImputeZeroes(innerFold$test_data,impFactor = innerFold$imp_factor)) 
      
      ## compute log ratios
      lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = innerFold$y_train,trainx))
      lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = innerFold$y_test,testx))
      
      cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train, 
                                          includeInfoGain = T, nfolds = 1, numRepeats = 1, 
                                          rankOrder = F)
      
    }))

    
    
    for(tar_Features in sets){
      
      suppressMessages(suppressWarnings({
      
        tar_dcvInner = targeted_dcvSelection(trainx = trainx,
                                        testx = testx,
                                        dcv = cc.dcv$lrs,lrs.train = lrs.train,lrs.test = lrs.test,
                                        y_label = innerFold$y_train,
                                        seed = sd1,
                                        ensemble = ensemble,
                                        y_test = innerFold$y_test,
                                        tarFeatures = tar_Features,
                                        ts.id = innerFold$test_ids, 
                                        max_sparsity = max_sparsity
        )
        
        perf = data.frame(Seed = sd1,Fold = f,tar_Features ,tar_dcvInner$Performance)
        inner_perf = rbind(inner_perf,perf)
      }))
      
      message(tar_Features)
      
  }
  
}

}

## aggregate results
inner_perf1 = inner_perf %>% 
  dplyr::group_by(Approach,tar_Features) %>% 
  summarise_all(.funs = mean)
ggplot(inner_perf1,aes(tar_Features,AUC,col = Approach))+
  geom_point()+
  geom_line()

inner_perf2 = inner_perf %>% 
  dplyr::group_by(tar_Features) %>% 
  summarise_all(.funs = mean)
ggplot(inner_perf2,aes(tar_Features,AUC))+
  geom_point()+
  geom_line()








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
tar_Features = 67
tar_Features = inner_perf2$tar_Features[which.max(inner_perf2$AUC)]

  
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

perf = data.frame(Scenario = "dcv_targeted",Dataset = f_name,tar_Features,
                  Seed = sd,Fold = f,Partition = "Test",tar_dcv$Performance)
benchmark = rbind(benchmark,perf)

tar_dcv$weighted_features$train

