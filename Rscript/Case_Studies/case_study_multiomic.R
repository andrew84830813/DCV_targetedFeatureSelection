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
# clus <- parallel::makeCluster(10)
# doParallel::registerDoParallel(clus)




# Read External Args ---------------------------------------------------------------

args = c(2,5,0)
args = commandArgs(trailingOnly = TRUE)
sd = as.numeric(args[1]) # random seed selection
f = as.numeric(args[2])
permute_labels = as.logical(as.numeric(args[3])) #should be 0(False) or 1(True)


## set random seed
set.seed(sd)

f_name = paste0("ihmp-multiomics_seed",sd,"_fold",f)

# Partition Data ----------------------------------------------------------

set.seed(sd)
load("Output/ihmp/folds.Rda")
fld = folds_seeds[[sd]]$folds
k_fold = max(fld)


## load multiomic data
load("Output/ihmp/ihmp_multiomic_data.Rda")
data_types = names(multiomic)  
multiomic.train = data.frame()
multiomic.test = data.frame()
benchmark = data.frame()



  for(dt in data_types){
    
    df = multiomic[[dt]]
    
    ## Partition Data
    allData = lodo_partition(df,dataset_labels = fld,sd)
    
    
    ## Extract Test/Train Split
    ttData = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                              fold = f,
                                              maxSparisty = .9,
                                              extractTelAbunance = F)
    ##get train test partitions
    train.data = ttData$train_Data
    test.data = ttData$test_data
    ## Compute Total Parts
    number_parts = ncol(train.data);number_parts
    nsamps = nrow(train.data)
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
    #ensemble = c("ranger","xgbTree","gbm")
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
    # ggplot(inner_perf1,aes(tar_Features,AUC,col = Approach))+
    #   geom_point()+
    #   geom_line()
    inner_perf2 = inner_perf %>% 
      dplyr::group_by(tar_Features) %>% 
      summarise_all(.funs = mean)
    # ggplot(inner_perf2,aes(tar_Features,AUC))+
    #   geom_point()+
    #   geom_line()
    
    
    ## Train final Model
    trainx = data.frame(fastImputeZeroes(ttData$train_Data,impFactor = ttData$imp_factor))
    testx = data.frame(fastImputeZeroes(ttData$test_data,impFactor = ttData$imp_factor))
    ## compute log ratios
    lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = ttData$y_train,trainx))
    lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = ttData$y_test,testx))
    ## Comput DCV
    cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train, 
                                        includeInfoGain = T, nfolds = 1, numRepeats = 1, 
                                        rankOrder = F)
    ## Define Target feature from nested cross validation
    tar_Features = inner_perf2$tar_Features[which.max(inner_perf2$AUC)]
    ## Run Targetd Feature Selection
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
    
    ## Compute Multiclass AUC statistics
    apps = names(tar_dcv$probMatrices)
    mulit_perf = data.frame()
    for(ap  in apps){
      pr =tar_dcv$probMatrices[[ap]]
      mroc = pROC::multiclass.roc( pr$Status, pr[,classes]);mroc
      c= combinat::combn2(classes)
      ## multiclass pairwise AUC
      pwAUC = data.frame()
      for(jid in 1:nrow(c) ){
        pw = as.character(c[jid,])
        prn = pr
        prn = prn[prn$Status%in%pw,]
        mc = unique(prn$Status)
        if(sum(mc%in%c("nonIBD"))>1){
          controls = which(mc%in%c("nonIBD"))
          cases = which(!mc%in%c("nonIBD"))
          mc = c(mc[controls],mc[cases])
        }else{
          mc = as.character(unique(prn$Status))
        }
        
        optThreshold = pROC::roc(prn$Status,prn[,mc[1]],levels = mc)
        thrhold = pROC::coords(optThreshold,x="best",
                               ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                               transpose = FALSE)
        optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
        kk = which.max(optMCC)
        oo = data.frame(ID = jid,Comp = paste0(pw[1],"_",pw[2]),mcc = optMCC[kk],thrhold[kk,],
                        AUC = as.numeric(pROC::auc(optThreshold)))
        pwAUC = rbind(pwAUC,oo)
      }
      colnames(pwAUC) = paste0("pw_",colnames(pwAUC))
      imp_ids = c(3,9:12,14)
      pairwise = data.frame(t(colMeans(pwAUC[,imp_ids],na.rm = T)))
      colnames(pairwise) = str_split(colnames(pairwise),pattern = "_",simplify = T)[,2]
      for(ind in imp_ids){
        nm = str_split(colnames(pwAUC)[ind],pattern = "_",simplify = T)[,2]
        thint = spread(pwAUC[,c(2,ind)],"pw_Comp",colnames(pwAUC)[ind])
        colnames(thint) = paste0(colnames(thint),paste0(".",nm))      
        pairwise =  cbind(pairwise,thint)
      } 
      pairwise$AUC = as.numeric(pROC::auc(mroc))
      ph = data.frame(Approach = ap,pairwise)
      mulit_perf = rbind(mulit_perf,ph)
    }
    c1 = which(colnames(mulit_perf)%in%c("UC_CD.AUC","CD_UC.AUC"))
    c2 = which(colnames(mulit_perf)%in%c("nonIBD_CD.AUC","CD_nonIBD.AUC"))
    c3 = which(colnames(mulit_perf)%in%c("nonIBD_UC.AUC","UC_nonIBD.AUC"))
    mulit_perf = data.frame(Approach = mulit_perf$Approach,MCC = mulit_perf$mcc,AUC = mulit_perf$AUC,
                            AUC.UC_CD = mulit_perf[,c1],
                            AUC.nonIBD_UC = mulit_perf[,c3],
                            AUC.nonIBD_CD = mulit_perf[,c2])
   
    ## Write Results
    perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),data_type = dt,
                      Dataset = f_name,Seed = sd,Fold = f,
                      number_parts = tar_dcv$Performance$number_parts,
                      number_ratios = tar_dcv$Performance$number_ratios ,
                      comp_time = tar_dcv$Performance$comp_time,
                      base_dims = ncol(train.data),mulit_perf
    )
    
    benchmark = rbind(benchmark,perf)
    
    ## Store data type specif features - train 
    merdat = data.frame(ttData$train_ids[,-3],outer_fold = f,seed = sd,tar_dcv$rfe_features$train)
    merdat = gather(merdat,"log_ratio","value",5:ncol(merdat))
    merdat$log_ratio = paste0(dt,"._.",merdat$log_ratio)
    multiomic.train = rbind(multiomic.train,merdat)
    ## Store data type specif features - test 
    merdat = data.frame(ttData$test_ids[,-3],outer_fold = f,seed = sd,tar_dcv$rfe_features$test)
    merdat = gather(merdat,"log_ratio","value",5:ncol(merdat))
    merdat$log_ratio = paste0(dt,"._.",merdat$log_ratio)
    multiomic.test = rbind(multiomic.test,merdat)
    
    
  }
  
  
  ## Merge All features across data types - Spread data
  multiomic.train1 = spread(multiomic.train,key = "log_ratio",value = "value")
  multiomic.test1 = spread(multiomic.test,key = "log_ratio",value = "value")
  ## retrieve test and train data
  train.data = multiomic.train1[,-4:-1]
  test.data = multiomic.test1[,-4:-1]
  y_label = multiomic.train1[,2]
  y_test = multiomic.test1[,2]
  
  
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
  pmat = data.frame(p[,,1])
  pr = data.frame(Status = y_test,pmat)
  mroc = pROC::multiclass.roc( pr$Status, pr[,classes]);mroc
  c= combinat::combn2(classes)
  ## multiclass pairwise AUC
  pwAUC = data.frame()
  for(jid in 1:nrow(c) ){
    pw = as.character(c[jid,])
    prn = pr
    prn = prn[prn$Status%in%pw,]
    mc = unique(prn$Status)
    if(sum(mc%in%c("nonIBD"))>1){
      controls = which(mc%in%c("nonIBD"))
      cases = which(!mc%in%c("nonIBD"))
      mc = c(mc[controls],mc[cases])
    }else{
      mc = as.character(unique(prn$Status))
    }
    
    optThreshold = pROC::roc(prn$Status,prn[,mc[1]],levels = mc)
    thrhold = pROC::coords(optThreshold,x="best",
                           ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                           transpose = FALSE)
    optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
    kk = which.max(optMCC)
    oo = data.frame(ID = jid,Comp = paste0(pw[1],"_",pw[2]),mcc = optMCC[kk],thrhold[kk,],
                    AUC = as.numeric(pROC::auc(optThreshold)))
    pwAUC = rbind(pwAUC,oo)
  }
  colnames(pwAUC) = paste0("pw_",colnames(pwAUC))
  imp_ids = c(3,9:12,14)
  pairwise = data.frame(t(colMeans(pwAUC[,imp_ids],na.rm = T)))
  colnames(pairwise) = str_split(colnames(pairwise),pattern = "_",simplify = T)[,2]
  for(ind in imp_ids){
    nm = str_split(colnames(pwAUC)[ind],pattern = "_",simplify = T)[,2]
    thint = spread(pwAUC[,c(2,ind)],"pw_Comp",colnames(pwAUC)[ind])
    colnames(thint) = paste0(colnames(thint),paste0(".",nm))      
    pairwise =  cbind(pairwise,thint)
  } 
  pairwise$AUC = as.numeric(pROC::auc(mroc))
  ph = data.frame(Approach = "penRegression",pairwise)
  c1 = which(colnames(ph)%in%c("UC_CD.AUC","CD_UC.AUC"))
  c2 = which(colnames(ph)%in%c("nonIBD_CD.AUC","CD_nonIBD.AUC"))
  c3 = which(colnames(ph)%in%c("nonIBD_UC.AUC","UC_nonIBD.AUC"))
  ph = data.frame(Approach = ph$Approach,MCC = ph$mcc,AUC = ph$AUC,
                          AUC.UC_CD = ph[,c1],
                          AUC.nonIBD_UC = ph[,c3],
                          AUC.nonIBD_CD = ph[,c2])
  
  cn = feat.df$Ratio
  n_ratios = length(cn)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)
  
  perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                    data_type = "integrated-penRegression",
                    Dataset = f_name,Seed = sd,Fold = f,
                    number_parts = n_parts,
                    number_ratios = nrow(feat.df) ,
                    comp_time = NA,
                    base_dims = ncol(train.data),ph
  )
  
  benchmark = rbind(benchmark,perf)
  
  
  ## Use RFE To select features and train model
  compTime2 = system.time({
    pp = rfeSelection.ByMetric(train_ratio = train.data,
                               test_ratio = test.data,
                               ytrain =y_label,
                               ntrees = 2000,
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
                      y_test = y_test,
                      testIDs = multiomic.test1[,1:4],
                      models = ensemble) 
  
  ## Compute Performance
  
  pmat = ph$predictionMatrix
  pmat = pmat %>% 
    group_by(ID,Status) %>% 
    dplyr::select(-model) %>% 
    summarise_all(.funs = mean)
  pmat = data.frame(pmat)
  mroc = pROC::multiclass.roc(pmat$Status,pmat[,classes]);mroc
  
  
  ## Compute Number of Part
  cn = colnames(train_data2)
  n_ratios = length(cn)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)
  
 
    pr =pmat
    mroc = pROC::multiclass.roc( pr$Status, pr[,classes]);mroc
    c= combinat::combn2(classes)
    ## multiclass pairwise AUC
    pwAUC = data.frame()
    for(jid in 1:nrow(c) ){
      pw = as.character(c[jid,])
      prn = pr
      prn = prn[prn$Status%in%pw,]
      mc = unique(prn$Status)
      if(sum(mc%in%c("nonIBD"))>1){
        controls = which(mc%in%c("nonIBD"))
        cases = which(!mc%in%c("nonIBD"))
        mc = c(mc[controls],mc[cases])
      }else{
        mc = as.character(unique(prn$Status))
      }
      
      optThreshold = pROC::roc(prn$Status,prn[,mc[1]],levels = mc)
      thrhold = pROC::coords(optThreshold,x="best",
                             ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                             transpose = FALSE)
      optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
      kk = which.max(optMCC)
      oo = data.frame(ID = jid,Comp = paste0(pw[1],"_",pw[2]),mcc = optMCC[kk],thrhold[kk,],
                      AUC = as.numeric(pROC::auc(optThreshold)))
      pwAUC = rbind(pwAUC,oo)
    }
    colnames(pwAUC) = paste0("pw_",colnames(pwAUC))
    imp_ids = c(3,9:12,14)
    pairwise = data.frame(t(colMeans(pwAUC[,imp_ids],na.rm = T)))
    colnames(pairwise) = str_split(colnames(pairwise),pattern = "_",simplify = T)[,2]
    for(ind in imp_ids){
      nm = str_split(colnames(pwAUC)[ind],pattern = "_",simplify = T)[,2]
      thint = spread(pwAUC[,c(2,ind)],"pw_Comp",colnames(pwAUC)[ind])
      colnames(thint) = paste0(colnames(thint),paste0(".",nm))      
      pairwise =  cbind(pairwise,thint)
    } 
    pairwise$AUC = as.numeric(pROC::auc(mroc))
    ph = data.frame(Approach = "rfe",pairwise)
    c1 = which(colnames(ph)%in%c("UC_CD.AUC","CD_UC.AUC"))
    c2 = which(colnames(ph)%in%c("nonIBD_CD.AUC","CD_nonIBD.AUC"))
    c3 = which(colnames(ph)%in%c("nonIBD_UC.AUC","UC_nonIBD.AUC"))
    ph = data.frame(Approach = ph$Approach,MCC = ph$mcc,AUC = ph$AUC,
                    AUC.UC_CD = ph[,c1],
                    AUC.nonIBD_UC = ph[,c3],
                    AUC.nonIBD_CD = ph[,c2])
    
  
  perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                    data_type = "integrated-rfe",
                    Dataset = f_name,Seed = sd,Fold = f,
                    number_parts = n_parts,
                    number_ratios = ncol(train_data2) ,
                    comp_time = NA,
                    base_dims = ncol(train.data),ph
  )
  
  benchmark = rbind(benchmark,perf)
  
  
  
  write_csv(x = benchmark,file = paste0("Results/",f_name,".csv"))
  
