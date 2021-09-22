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
library(selbal) # selbal



#setwd("/nas/longleaf/home/andrew84/rProjects/DiCoVarFS_project")


# Load Helper Functions  ---------------------------------------------------------------
fnames = dir("Helper_Functions/")
for(f in fnames){
  source(paste0("Helper_Functions/",f))
}



# Setup Cluster ---------------------------------------------------------------

## Detect Cores
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)




# Read External Args ---------------------------------------------------------------

args = c(2,5,0)
args = commandArgs(trailingOnly = TRUE)
sd = as.numeric(args[1]) # random seed selection
f = as.numeric(args[2])
permute_labels = as.logical(as.numeric(args[3])) #should be 0(False) or 1(True)


## set random seed
set.seed(sd)
###-------------------------------------------------------*
#### NEC ####
###-------------------------------------------------------*
metaData = data.frame(readr::read_tsv(file = "Data/NICUNEC.WGS.sample_details.tsv"))
otuData = data.frame(readr::read_tsv(file = "Data/NICUNEC.WGS.taxon_abundance.tsv"))


## Filter NEC samples for NEC  in the that occurs a least a week from now
metaData1 = metaData %>% 
  filter(Days.of.period.NEC.diagnosed<=-7 ) %>% 
  filter(Age.at.sample.collection..days.>0) 

## Control Samples
metaData2 = metaData %>% 
  filter(is.na(Days.of.period.NEC.diagnosed)) %>% 
  filter(Age.at.sample.collection..days. < max(metaData1$Age.at.sample.collection..days.)) 
hist(metaData2$Age.at.sample.collection..days.)
hist(metaData1$Age.at.sample.collection..days.)


## Merge Samples
set.seed(08272008)
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





benchmark = data.frame()
permute_labels = F
f_name= "case_study_NEC"

for(sd in 1:10){
  set.seed(sd)
  
  ## Partition Data
  k_fold = 5
  ##strat by sample and label
  overll_folds = caret::createFolds(samples$Necrotizing.enterocolitis,k =k_fold,list = F)
  samples1 = samples
  samples1$fold = overll_folds
  samples1 = left_join(md1,samples1)
  
  allData = lodo_partition(df,dataset_labels = samples1$fold,sd)
  
  message("\n\n``````` Start Seed:  ",sd,"````````````\n\n")
  
  #  Perform 2-fold cross validation -------------------------------
  
  for(f in 1:k_fold){
    
    ## Extract Test/Train Splilt
    ttData = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                              fold = f,
                                              maxSparisty = .9,
                                              extractTelAbunance = F)
    ##get train test paritions
    train.data = ttData$train_Data
    test.data = ttData$test_data
    ## Compute Total Parts
    number_parts = ncol(train.data);number_parts
    nsamps = nrow(train.data)
    table(ttData$y_test)
    classes = as.character(unique(ttData$y_train))
    
    ##get metadata
    train.md = data.frame(X1 = rownames(train.data))
    train.md = left_join(train.md,md1)
    test.md = data.frame(X1 = rownames(test.data))
    test.md = left_join(test.md,md1)
    y_label = ttData$y_train

    
   ## Pre-Process Metadata
    mm.char = train.md %>% 
      select(Sex,
             Chorioamnionitis,
             Maternal.antepartum.antibiotics.administered,
             Baby.s.delivery.mode,
             Born.stunted,
             Host.diet)
    mm.char.test = subset(test.md,select = colnames(mm.char))
    ##-----------------------------
    ## Add covairate
    ## cat features
    dummy <- dummyVars(" ~ .", data=rbind(mm.char,mm.char.test))
    newdata <- data.frame(predict(dummy, newdata = rbind(mm.char,mm.char.test)))
    newdata.train = newdata[1:nrow(mm.char),]
    newdata.test =  newdata[-1:-nrow(mm.char),]
    
    
    ## cont features
    train_metadata = data.frame(Age = train.md$Age.at.sample.collection..days.,
                               gesAge = train.md$Gestational.age.at.birth..days.,
                               host_weight = train.md$Host.weight..g.,
                               newdata.train
    )
    
    
    test.metadata = data.frame(Age = test.md$Age.at.sample.collection..days.,
                               gesAge = test.md$Gestational.age.at.birth..days.,
                               host_weight = test.md$Host.weight..g.,
                               newdata.test
    )
    
    
    
  
    
    
    # RFE - Metadata --------------------------------------------------------------------
    ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")
    suppressMessages(suppressWarnings({
  
      compTime2 = system.time({
        pp = rfeSelection.ByMetric(train_ratio = train_metadata,
                                   test_ratio = test.metadata,
                                   ytrain =y_label,
                                   ntrees = 750,
                                   sets = 10,
                                   impMeasure = "impurity_corrected",
                                   kfold = 5,
                                   minPercentFeatReturn = .3)
      })
      train_data.meta = pp$reducedTrainRatios
      test_data.meta = pp$reducedTestRatio
      message("number of features = ",ncol(train_data.meta))
      keep.meta = colnames(train_data.meta)
      
      
      ## Train Model
      ph = trainML_Models(trainLRs = train_data.meta,
                          testLRs = test_data.meta,
                          ytrain = y_label,
                          y_test = ttData$y_test,
                          testIDs = ttData$test_ids,
                          models = ensemble) 
      
      ## Compute Performance
      geo.mean = function(x){
        exp(mean(log(x)))
      }
      pmat = ph$predictionMatrix
      pmat = pmat %>% 
        group_by(ID,Status) %>% 
        dplyr::select(-model) %>% 
        summarise_all(.funs = mean)
      pmat = data.frame(pmat)
      mroc = pROC::multiclass.roc(pmat$Status,pmat[,classes]);mroc
      
      
      ## Compute Number of Part
      cn = colnames(train_data.meta)
      n_ratios = length(cn)
      uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
      n_parts  = dplyr::n_distinct(uniqueParts)
      
      ## Save Performance
      perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                        Dataset = f_name,Seed = sd,Fold = f,Approach = "meta_alone",AUC = as.numeric(pROC::auc(mroc)),
                        number_parts = n_parts,number_ratios = ncol(train_data.meta) ,comp_time = NA,
                        base_dims = ncol(train.data))
      benchmark = rbind(benchmark,perf)
    }))
    
    
    # RFE - Mbiome --------------------------------------------------------------------
    
    perc_totalParts2Keep = .75
    num_sets = 5
    
    base_dims = ncol(ttData$train_Data)
    max_parts = round(perc_totalParts2Keep*base_dims)
    sets = round(seq(1,max_parts,length.out = num_sets))[-1]
    
    
    inner_perf = data.frame()
    
    
    ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")
    #ensemble = c("ranger","xgbTree","xgbLinear")
    max_sparsity = .9
    
    
    ## Tune target features
    for(sd1 in 1:1){
      set.seed(sd1)
      k_fold = 2
      overll_folds = caret::createFolds(ttData$y_train,k = k_fold,list = F)
      innerfold_data = lodo_partition(data.frame(Status = ttData$y_train,ttData$train_Data),
                                      dataset_labels = overll_folds,
                                      sd1)
      
      
      ## Get within fold cross validated performance 
      
      for(ff in 1:k_fold){
        
        ## Partition inner fold
        innerFold = DiCoVarML::extractTrainTestSplit(foldDataList = innerfold_data,
                                                     fold = ff,
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
            
            perf = data.frame(Seed = sd1,Fold = ff,tar_Features ,tar_dcvInner$Performance)
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
    
    perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                      Dataset = f_name,Seed = sd,Fold = f,Approach = tar_dcv$Performance$Approach,AUC = as.numeric(tar_dcv$Performance$AUC),
                      number_parts = tar_dcv$Performance$number_parts,number_ratios = tar_dcv$Performance$number_ratios ,
                      comp_time = tar_dcv$Performance$comp_time,
                      base_dims = ncol(train.data)
    )
    benchmark = rbind(benchmark,perf)
    
    
    
    
    # RFE - Mbiome + Meta --------------------------------------------------------------------
    suppressMessages(suppressWarnings({
     
      
      train_data2 = cbind(tar_dcv$rfe_features$train,train_data.meta)
      test_data2 = cbind(tar_dcv$rfe_features$test,test_data.meta)
      
      
      compTime2 = system.time({
        pp = rfeSelection.ByMetric(train_ratio = train_data2,
                                   test_ratio = test_data2,
                                   ytrain =y_label,
                                   ntrees = 750,
                                   sets = 10,
                                   impMeasure = "impurity_corrected",
                                   kfold = 5,
                                   minPercentFeatReturn = .3)
      })
      train_data2 = pp$reducedTrainRatios
      test_data2 = pp$reducedTestRatio
      message("number of features = ",ncol(train_data2))
      
      ## Train Model
      ph = trainML_Models(trainLRs = train_data2,
                          testLRs = test_data2,
                          ytrain = y_label,
                          y_test = ttData$y_test,
                          testIDs = ttData$test_ids,
                          models = ensemble) 
      
      ## Compute Performance
      geo.mean = function(x){
        exp(mean(log(x)))
      }
      pmat = ph$predictionMatrix
      pmat = pmat %>% 
        group_by(ID,Status) %>% 
        dplyr::select(-model) %>% 
        summarise_all(.funs = mean)
      pmat = data.frame(pmat)
      #pmat[,4:ncol(pmat)] = clo((pmat[,4:ncol(pmat)]))
      mroc = pROC::multiclass.roc(pmat$Status,pmat[,classes]);mroc
      
      
      ## Compute Number of Part
      cn = colnames(train_data2)
      n_ratios = length(cn)
      uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
      n_parts  = dplyr::n_distinct(uniqueParts)
      
      ## Save Performance
      perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                        Dataset = f_name,Seed = sd,Fold = f,Approach = "mbiome + meta_data",AUC = as.numeric(pROC::auc(mroc)),
                        number_parts = n_parts,number_ratios = ncol(train_data2) ,comp_time = NA,
                        base_dims = ncol(train.data))
      benchmark = rbind(benchmark,perf)
    }))
    
    
    
    
   
  }
}

write_csv(benchmark,file = "Results/case_studyNEC_results.csv")


f <- function(x) {
  r <- quantile(x, probs = c(0.10, 0.25, 0.5, 0.75, 0.90))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

 res = benchmark %>% 
  group_by(Approach,Dataset) %>% 
  summarise_all(mean)
res = data.frame(res)

res = benchmark %>% 
  group_by(Approach,Dataset,Seed) %>% 
  summarise_all(mean)


res1 = benchmark %>% 
  group_by(Approach,Dataset) %>% 
  summarise_all(mean)

ggplot(res,aes(Approach,AUC))+
  stat_summary(fun.y  = mean,geom = "point")+
  stat_summary(fun.data = mean_se,geom = "errorbar")
  


## kruskal test
kw = spread(benchmark[,1:6],"Approach","AUC")
kruskal.test(x = res$AUC,g = res$Approach)

wilcox.test(x = kw$`mbiome + meta_data`,y = kw$`DCV-rfRFE`,paired = T,alternative = "two.sided")
wilcox.test(x = kw$`DCV-ridgeEnsemble`,y = kw$`Coda-LASSO`,paired = T,alternative = "two.sided")


