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

args = c(2,3,0)
args = commandArgs(trailingOnly = TRUE)
sd = as.numeric(args[1]) # random seed selection
f = as.numeric(args[2]) # random seed selection
permute_labels = as.logical(as.numeric(args[3])) #should be 0(False) or 1(True)


## set random seed
set.seed(sd)

## Name Iteration
f_name = paste0("NEC_seed",sd,"_permute",permute_labels,"_fold",f)

## read data
load("Output/microbiomeDB_dataAtLeast7DaysBeforeNEC.Rda")
df = exp$mbiome
md1 = exp$metadata

## Get fold matrix
folds_matrix = data.frame(read_csv(file = "Output/commonFolds_NEC_crossValidation.csv"))


## Partition data
allData = lodo_partition(df,dataset_labels = folds_matrix[,sd],sd)

##Define Out Data frame
benchmark = data.frame()


message("\n\n``````` Start Seed:  ",sd,"````````````\n\n","fold",f)

## Compute Performance for the f - th fold

  ## Extract Test/Train Split
  ttData = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                            fold = f,
                                            maxSparisty = .9,
                                            permLabels = permute_labels,
                                            extractTelAbunance = F)
  ##get train test partitions
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
  message("Compute Performance - Metadata Alone RFE")
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
                      Dataset = f_name,Seed = sd,Fold = f,Approach = "meta_aloneRFE",AUC = as.numeric(pROC::auc(mroc)),
                      number_parts = n_parts,number_ratios = ncol(train_data.meta) ,comp_time = NA,
                      base_dims = ncol(train.data))
    benchmark = rbind(benchmark,perf)
  }))
  
  
  # GLM -  Meta --------------------------------------------------------------------
  message("Compute Performance - Metadata Alone GLM")
  suppressMessages(suppressWarnings({
    
    
    ## retrieve test and train data
    train.data = cbind(train_metadata)
    test.data = cbind(test.metadata)
    y_label = ttData$y_train
    y_test = ttData$y_test
    
    
    ## Apply Penalized Regression
    ## Tune Alpha
    type_family = if_else(length(classes)>2,"multinomial","binomial")
    infld = 2
    flds = caret::createFolds(y = y_label,k = infld,list = F)
    compTime2 = system.time({
      aseq = seq(1e-3,1,length.out = 10)
      min_dev =  foreach(a = aseq,.combine = rbind)%dopar%{
        
        aucc = c()
        for(f in 1:infld){
          bool  = flds==f
          compTime2 = system.time({
            cv.clrlasso <- glmnet::cv.glmnet(as.matrix(train.data[bool,]),y_label[bool], 
                                             standardize=T, alpha=a,family=type_family)
          })
          
          ## make predictions
          p = predict(cv.clrlasso, newx = as.matrix(train.data[!bool,]), s = "lambda.min",type = "response")
          if(type_family=="binomial"){
            mroc = pROC::roc(y_label[!bool],p)
            mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
          }else{
            ## multiclass
            mroc = pROC::multiclass.roc(y_label[!bool],p[,,1])
            mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
          }
          aucc = c(aucc,as.numeric(mroc.dcvlasso))
        }
        data.frame(a,auc = mean(aucc))
      }
    })
    min_dev = min_dev %>% 
      arrange(desc(auc))
    
    ## Train GLM
    compTime2 = system.time({
      cv.clrlasso <- glmnet::cv.glmnet(as.matrix(train.data),y_label, standardize=T, alpha = min_dev$a[1],family=type_family)
    })
    if(type_family=="binomial"){
      features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
      features = features[-1,]
      features = features[abs(features)>0]
      length(features)
      c = as.matrix(coef(cv.clrlasso, s = "lambda.min"))[-1,]
      train_data.metaGLM = subset(train.data,select = names(features))
      test_data.metaGLM = subset(test.data,select = names(features))
    }else{
      features = as.matrix(stats::coef(cv.clrlasso, s = "lambda.min"))
      feat.df = data.frame()
      for(o in 1:length(features)){
        ph = as.matrix(features[[o]])
        feat = ph[-1,]
        keep = feat[abs(feat)>0]
        feat.df = rbind(feat.df,data.frame(Ratio = names(keep),coef = as.numeric(keep)))
      }
      feat.df =feat.df %>% 
        group_by(Ratio) %>% 
        summarise(coef = sum(coef)) %>% 
        filter(coef!=0)
      train_data.metaGLM = subset(train.data,select = feat.df$Ratio)
      test_data.metaGLM = subset(test.data,select = feat.df$Ratio)
    }
    
    ## make predictions
    p = predict(cv.clrlasso, newx = as.matrix(test.data), s = "lambda.min",type = "response")
    if(type_family=="binomial"){
      mroc = pROC::roc(ttData$y_test,p)
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
    }else{
      ## multiclass
      mroc = pROC::multiclass.roc(ttData$y_test,p[,,1])
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
    }
    
    ## Compute Number of Part
    cn =names(c[abs(c)>0])
    n_ratios = length(cn)
    uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
    n_parts  = dplyr::n_distinct(uniqueParts)
    
    ## Save Performance
    perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                      Dataset = f_name,Seed = sd,Fold = f,Approach = "meta_dataGLM",
                      AUC = as.numeric(pROC::auc(mroc)),
                      number_parts = n_parts,number_ratios = ncol(train.data) ,comp_time = NA,
                      base_dims = ncol(train.data))
    benchmark = rbind(benchmark,perf)
    
    
  }))
  
  
  # RFE - Mbiome --------------------------------------------------------------------
  message("Compute Performance - Mbiome Target DCV")
  suppressMessages(suppressWarnings({
    
    perc_totalParts2Keep = .75
    num_sets = 5
    
    base_dims = ncol(ttData$train_Data)
    max_parts = round(perc_totalParts2Keep*base_dims)
    sets = round(seq(1,max_parts,length.out = num_sets))[-1]
    
    
    inner_perf = data.frame()
    
    
    ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")
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
                      Dataset = f_name,Seed = sd,Fold = f,Approach = paste0("mbiome_",tar_dcv$Performance$Approach),AUC = as.numeric(tar_dcv$Performance$AUC),
                      number_parts = tar_dcv$Performance$number_parts,number_ratios = tar_dcv$Performance$number_ratios ,
                      comp_time = tar_dcv$Performance$comp_time,
                      base_dims = ncol(train.data)
    )
    benchmark = rbind(benchmark,perf)
  
  }))
  
  
  # RFE - Mbiome + Meta --------------------------------------------------------------------
  message("Compute Performance - Mbiome + Metadata RFE")
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
                      Dataset = f_name,Seed = sd,Fold = f,Approach = "mbiome + meta_dataRFE",AUC = as.numeric(pROC::auc(mroc)),
                      number_parts = n_parts,number_ratios = ncol(train_data2) ,comp_time = NA,
                      base_dims = ncol(train.data))
    benchmark = rbind(benchmark,perf)
  }))
  
  
  
  # GLM - Mbiome + Meta --------------------------------------------------------------------
  message("Compute Performance - Mbiome + Metadata GLM")
  suppressMessages(suppressWarnings({
    
   
    ## retrieve test and train data
    cc = tar_dcv$ridge_coefficients
    train.data = cbind(sweep(tar_dcv$weighted_features$train,MARGIN = 2,STATS = as.numeric(cc),FUN = "*")  ,train_data.metaGLM)
    test.data = cbind(sweep(tar_dcv$weighted_features$test,MARGIN = 2,STATS = as.numeric(cc),FUN = "*"),test_data.metaGLM)
    y_label = ttData$y_train
    y_test = ttData$y_test
    
    
    ## Apply Penalized Regression
    ## Tune Alpha
    type_family = if_else(length(classes)>2,"multinomial","binomial")
    infld = 2
    flds = caret::createFolds(y = y_label,k = infld,list = F)
    compTime2 = system.time({
      aseq = seq(1e-3,1,length.out = 10)
      min_dev =  foreach(a = aseq,.combine = rbind)%dopar%{
        
        aucc = c()
        for(f in 1:infld){
          bool  = flds==f
          compTime2 = system.time({
            cv.clrlasso <- glmnet::cv.glmnet(as.matrix(train.data[bool,]),y_label[bool], 
                                             standardize=T, alpha=a,family=type_family)
          })
          
          ## make predictions
          p = predict(cv.clrlasso, newx = as.matrix(train.data[!bool,]), s = "lambda.min",type = "response")
          if(type_family=="binomial"){
            mroc = pROC::roc(y_label[!bool],p)
            mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
          }else{
            ## multiclass
            mroc = pROC::multiclass.roc(y_label[!bool],p[,,1])
            mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
          }
          aucc = c(aucc,as.numeric(mroc.dcvlasso))
        }
        data.frame(a,auc = mean(aucc))
      }
    })
    min_dev = min_dev %>% 
      arrange(desc(auc))
    
    ## Train GLM
    compTime2 = system.time({
      cv.clrlasso <- glmnet::cv.glmnet(as.matrix(train.data),y_label, standardize=T, alpha = min_dev$a[1],family=type_family)
    })
    if(type_family=="binomial"){
      features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
      features = features[-1,]
      features = features[abs(features)>0]
      length(features)
      c = as.matrix(coef(cv.clrlasso, s = "lambda.min"))[-1,]
    }else{
      features = as.matrix(stats::coef(cv.clrlasso, s = "lambda.min"))
      feat.df = data.frame()
      for(o in 1:length(features)){
        ph = as.matrix(features[[o]])
        feat = ph[-1,]
        keep = feat[abs(feat)>0]
        feat.df = rbind(feat.df,data.frame(Ratio = names(keep),coef = as.numeric(keep)))
      }
      feat.df =feat.df %>% 
        group_by(Ratio) %>% 
        summarise(coef = sum(coef)) %>% 
        filter(coef!=0)
    }
    
    ## make predictions
    p = predict(cv.clrlasso, newx = as.matrix(test.data), s = "lambda.min",type = "response")
    if(type_family=="binomial"){
      mroc = pROC::roc(ttData$y_test,p)
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
    }else{
      ## multiclass
      mroc = pROC::multiclass.roc(ttData$y_test,p[,,1])
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
    }
    
    ## Compute Number of Part
    cn =names(c[abs(c)>0])
    n_ratios = length(cn)
    uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
    n_parts  = dplyr::n_distinct(uniqueParts)
    
    ## Save Performance
    perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                      Dataset = f_name,Seed = sd,Fold = f,Approach = "mbiome + meta_dataGLM",
                      AUC = as.numeric(pROC::auc(mroc)),
                      number_parts = n_parts,number_ratios = ncol(train.data) ,comp_time = NA,
                      base_dims = ncol(train.data))
    benchmark = rbind(benchmark,perf)
    
    
  }))
   


## Write Performance Estimates
write_csv(x = benchmark,file = paste0("Results/",f_name,".csv"))
    


