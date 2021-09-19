
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


setwd("/nas/longleaf/home/andrew84/rProjects/DiCoVarFS_project/")

# Load Helper Functions  ---------------------------------------------------------------
fnames = dir("Helper_Functions/")
for(f in fnames){
  source(paste0("Helper_Functions/",f))
}


# Setup Cluster ---------------------------------------------------------------

## Detect Cores
#number_cores = parallel::detectCores()
# clus <- parallel::makeCluster(10)
# doParallel::registerDoParallel(clus)


# Read External Args ---------------------------------------------------------------

args = c(6,3,6)
args = commandArgs(trailingOnly = TRUE)
sd = as.numeric(args[1]) # random seed selection
f = as.numeric(args[2])
dataset = as.numeric(args[3])



# Load Data ---------------------------------------------------------------

## Read Data
switch (dataset,
        {load("Output/ihmp/proteomics_data.Rda");fname = "proteomics"},
        {load("Output/ihmp/metabolomics_C18-neg_data.Rda");fname = "C18-neg"},
        {load("Output/ihmp/metabolomics_C8-pos_data.Rda");fname = "C8-pos"},
        {load("Output/ihmp/metabolomics_HILIC-neg_data.Rda");fname = "HILIC-neg"},
        {load("Output/ihmp/metabolomics_HILIC-pos_data.Rda");fname = "HILIC-pos"},
        {load("Output/ihmp/mgxFunctions_data.Rda");fname = "mgxFunctions"},
        {load("Output/ihmp/metatranscriptomics_data.Rda");fname = "metatranscriptomics"},
        {load("Output/ihmp/mgxTaxonomy_data.Rda");fname = "mgxTaxonomy"}
)

## Which samples to keep 
load(file = "Output/ihmp/multiomic_keep.Rda")

## Get Data
df = obj$data
md = obj$md

##Select Key Samples
keep  =data.frame(External.ID = keep_samples$External.ID)
df = data.frame(External.ID = rownames(df),df)
df = left_join(keep,df)[,-1]
rownames(df) = keep$External.ID
md = left_join(keep,md)
## View Class Labels
table(df$Status)


## How many Taxa/Parts
message("Number of Taxa = ",ncol(df[,-1]))

## Degree of Sparsity
num_samples = n_distinct(md$Participant.ID)
subjectID = unique(md$Participant.ID)
totalParts = ncol(df[,-1])




# Partition Data ----------------------------------------------------------

set.seed(sd)
load("Output/ihmp/folds.Rda")
fld = folds_seeds[[sd]]$folds
## add folds and partition
df$fold = fld
df_ = df
df_ = data.frame(ID = rownames(df_), df_)
allData = list()
ph_list = list()
for (j in 1:dplyr::n_distinct(fld)) {
  xtrain = df %>% dplyr::filter(fold != j) %>% dplyr::select(-fold)
  xtrain_ID = df_ %>% dplyr::filter(fold != j) %>% dplyr::select(ID, 
                                                                 1, 2, fold)
  xtest = df %>% dplyr::filter(fold == j) %>% dplyr::select(-fold)
  xtest_ID = df_ %>% dplyr::filter(fold == j) %>% dplyr::select(ID, 
                                                                1, 2, fold)
  ph_list[[1]] = xtrain
  names(ph_list)[1] = "xtrain_combinedFolds"
  ph_list[[2]] = xtest
  names(ph_list)[2] = "xtest_kthFold"
  ph_list[[3]] = xtrain_ID
  names(ph_list)[3] = "xtrain_IDs"
  ph_list[[4]] = xtest_ID
  names(ph_list)[4] = "xtest_IDs"
  allData[[j]] = ph_list
  names(allData)[j] = paste("fold_", j, sep = "")
}
length(allData)



# Set Train Parms ----------------------------------------------------------

## Store Performance Data
perf.df = data.frame()
## Should Stacked Model Be Trained 
trainStacked = F
## Define Models/Ensemble
ensemble = c("ranger","pls","glmnet","svmRadial","gbm")
## Define Sparsity Threshold
max_sparsity = .9
## File info
fn = paste0("Results/benchmarkIHMP/Fold",f,"Seed",sd,"_Sparisty",max_sparsity*100,"_Dataset-",fname,"_resultsGLM.csv")
message(fn)



# Extract the kth-Fold Train/TEst Partition -------------------------------

ttData = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                          fold = f,
                                          maxSparisty = max_sparsity,
                                          extractTelAbunance = F)
## Compute Total Parts
number_parts = ncol(ttData$train_Data); message("Number of Parts: ",number_parts)
nsamps = nrow(ttData$train_Data)
table(ttData$y_test)
classes = as.character(unique(ttData$y_train))
if(length(classes)==2){
  controls = which(classes%in%c("Control","control","negative"))
  cases = which(!classes%in%c("Control","control","negative"))
  classes = c(classes[controls],classes[cases])
}


# GLM Selection  ---------------------------------------------------------------
selection_method = 2
methman = if_else(selection_method==1,"boruta","glm")


## Relative Abundance
message("GLM - Process Rel Abundance")
suppressMessages(suppressWarnings({
  compTime = system.time({
    model_data = DiCoVarML::relAbundanceFeatureSelection(train_data = ttData$train_Data,
                                                         y_train = ttData$y_train,
                                                         num_borutaRuns = 200,
                                                         test_data = ttData$test_data,
                                                         featureSelectionMethod = selection_method,
                                                         impute_factor = ttData$imputeFactor )
  })
  train = model_data$train_data
  test = model_data$test_data
  ## Train Model
  ph = trainML_Models(trainLRs =train,
                      testLRs =test,
                      ytrain = ttData$y_train,
                      y_test = ttData$y_test,
                      testIDs = ttData$test_ids,
                      models = ensemble ) 
  ## Extract Performance and Predictions
  p = ph$performance[,c(1:3,15)];colnames(p)[4] = "model"
  preds= ph$predictionMatrix
  preds = data.frame(Seed = sd,
                     Fold = f,
                     Sparsity = max_sparsity,
                     num_features = number_parts,
                     num_samples = nsamps,
                     num_parts = ncol(model_data$train_data) ,
                     num_ratios = NA,
                     method = paste0("relAbundance_",methman),
                     comp_time = compTime[1],
                     preds)
  preds = left_join(preds,p)
  mdl = unique(preds$model)
  for(m in mdl){
    pr = preds %>% 
      filter(model==m)
    mroc = pROC::multiclass.roc(pr$Status,pr[,classes])
    c= combinat::combn2(classes)
    ## multiclass pairwise AUC
    if(length(classes)>2){
      ## Compute Pairwise AUC
      pwAUC = data.frame()
      for(jid in 1:nrow(c) ){
        pw = as.character(c[jid,])
        prn = pr
        prn = prn[prn$Status%in%pw,]
        mc = unique(prn$Status)
        controls = which(mc%in%c("nonIBD","UC"))
        cases = which(!mc%in%c("nonIBD","UC"))
        mc = c(mc[controls],mc[cases])
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
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = ncol(model_data$train_data) ,
                      num_ratios = NA,
                      method =  paste0("relAbundance_",methman),
                      comp_time = compTime[1],
                      p[p$model==m,],pairwise)
      perf.df = rbind(perf.df,pp)
      
    }else{
      ## Compute Optimal Threshold and Metrics
      optThreshold = pROC::roc(pr$Status,pr[,classes[1]],levels = classes)
      thrhold = pROC::coords(optThreshold,x="best",
                             ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                             transpose = FALSE)
      optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
      kk = which.max(optMCC)
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = ncol(model_data$train_data) ,
                      num_ratios = NA,
                      method = paste0("relAbundance_",methman),
                      comp_time = compTime[1],mcc = optMCC[kk],thrhold[kk,],
                      p[p$model==m,],AUC = as.numeric(pROC::auc(mroc)))
      perf.df = rbind(perf.df,pp)
    }
  }
  
}))

## ALR
message("GLM - Process ALR")
suppressMessages(suppressWarnings({
  compTime = system.time({
    model_data = DiCoVarML::alrFeatureSelection(train_data = ttData$train_Data,
                                                y_train = ttData$y_train,num_borutaRuns = 200,
                                                test_data = ttData$test_data,
                                                featureSelectionMethod = selection_method,
                                                impute_factor = ttData$imputeFactor)
  })
  tr = model_data$train_data;colnames(tr) = paste0(colnames(tr),"___",colnames(ttData$train_Data)[1])
  ts = model_data$test_data;colnames(ts) = paste0(colnames(ts),"___",colnames(ttData$train_Data)[1])
  ## Train Model
  ph = trainML_Models(trainLRs = tr,
                      testLRs =ts,
                      ytrain = ttData$y_train,
                      y_test = ttData$y_test,testIDs = ttData$test_ids,models = ensemble ) 
  ## Extract Performance and Predictions
  cn = colnames(model_data$train_data)
  n_ratios = length(cn)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)
  p = ph$performance[,c(1:3,15)];colnames(p)[4] = "model"
  ## Predictions
  preds= ph$predictionMatrix
  preds = data.frame(Seed = sd,
                     Fold = f,
                     Sparsity = max_sparsity,
                     num_features = number_parts,
                     num_samples = nsamps,
                     num_parts = length(uniqueParts) ,
                     num_ratios = n_ratios,
                     method = paste0("alr_",methman),
                     comp_time = compTime[1],
                     preds)
  preds = left_join(preds,p)
  mdl = unique(preds$model)
  for(m in mdl){
    pr = preds %>% 
      filter(model==m)
    mroc = pROC::multiclass.roc(pr$Status,pr[,classes])
    
    
    c= combinat::combn2(classes)
    ## multiclass pairwise AUC
    if(length(classes)>2){
      ## Compute Pairwise AUC
      pwAUC = data.frame()
      for(jid in 1:nrow(c) ){
        pw = as.character(c[jid,])
        prn = pr
        prn = prn[prn$Status%in%pw,]
        mc = unique(prn$Status)
        controls = which(mc%in%c("nonIBD","UC"))
        cases = which(!mc%in%c("nonIBD","UC"))
        mc = c(mc[controls],mc[cases])
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
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = paste0("alr_",methman),
                      comp_time = compTime[1],
                      p[p$model==m,],pairwise)
      perf.df = rbind(perf.df,pp)
      
    }else{
      ## Compute Optimal Threshold and Metrics
      optThreshold = pROC::roc(pr$Status,pr[,classes[1]],levels = classes)
      thrhold = pROC::coords(optThreshold,x="best",
                             ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                             transpose = FALSE)
      optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
      kk = which.max(optMCC)
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = paste0("alr_",methman),
                      comp_time = compTime[1],mcc = optMCC[kk],thrhold[kk,],
                      p[p$model==m,],AUC = as.numeric(pROC::auc(mroc)))
      perf.df = rbind(perf.df,pp)
    }
    
    
    
    
  }
}))

## PLR
message("GLM - Process PLR")
suppressMessages(suppressWarnings({
  compTime = system.time({
    model_data = DiCoVarML::plrFeatureSelection(train_data = ttData$train_Data,
                                                y_train = ttData$y_train,num_borutaRuns = 200,
                                                test_data = ttData$test_data,
                                                featureSelectionMethod = selection_method,
                                                impute_factor = ttData$imputeFactor)
  })
  
  ## Train Model
  plr_featsGLM = model_data
  ph = trainML_Models(trainLRs = model_data$train_data,
                      testLRs =model_data$test_data,
                      ytrain = ttData$y_train,
                      y_test = ttData$y_test,
                      testIDs = ttData$test_ids,
                      models = ensemble ) 
  ## Extract Performance and Predictions
  cn = colnames(model_data$train_data)
  n_ratios = length(cn)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)
  p = ph$performance[,c(1:3,15)];colnames(p)[4] = "model"
  ## Predictions
  preds= ph$predictionMatrix
  preds = data.frame(Seed = sd,
                     Fold = f,
                     Sparsity = max_sparsity,
                     num_features = number_parts,
                     num_samples = nsamps,
                     num_parts = length(uniqueParts) ,
                     num_ratios = n_ratios,
                     method = paste0("plr_",methman),
                     comp_time = compTime[1],
                     preds)
  preds = left_join(preds,p)
  preds_plr = preds
  mdl = unique(preds$model)
  for(m in mdl){
    pr = preds %>% 
      filter(model==m)
    mroc = pROC::multiclass.roc(pr$Status,pr[,classes])
    
    
    c= combinat::combn2(classes)
    ## multiclass pairwise AUC
    if(length(classes)>2){
      ## Compute Pairwise AUC
      pwAUC = data.frame()
      for(jid in 1:nrow(c) ){
        pw = as.character(c[jid,])
        prn = pr
        prn = prn[prn$Status%in%pw,]
        mc = unique(prn$Status)
        controls = which(mc%in%c("nonIBD","UC"))
        cases = which(!mc%in%c("nonIBD","UC"))
        mc = c(mc[controls],mc[cases])
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
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = paste0("plr_",methman),
                      comp_time = compTime[1],
                      p[p$model==m,],pairwise)
      perf.df = rbind(perf.df,pp)
      
    }else{
      ## Compute Optimal Threshold and Metrics
      optThreshold = pROC::roc(pr$Status,pr[,classes[1]],levels = classes)
      thrhold = pROC::coords(optThreshold,x="best",
                             ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                             transpose = FALSE)
      optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
      kk = which.max(optMCC)
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = paste0("plr_",methman),
                      comp_time = compTime[1],mcc = optMCC[kk],thrhold[kk,],
                      p[p$model==m,],AUC = as.numeric(pROC::auc(mroc)))
      perf.df = rbind(perf.df,pp)
    }
    
    
  }
}))



# Boruta Selection  ---------------------------------------------------------------
selection_method = 1
methman = if_else(selection_method==1,"boruta","glm")

## Relative Abundance
message("Boruta - Process Rel Abundance")
suppressMessages(suppressWarnings({
  compTime = system.time({
    model_data = DiCoVarML::relAbundanceFeatureSelection(train_data = ttData$train_Data,
                                                         y_train = ttData$y_train,
                                                         num_borutaRuns = 200,
                                                         test_data = ttData$test_data,
                                                         featureSelectionMethod = selection_method,
                                                         impute_factor = ttData$imputeFactor )
  })
  train = model_data$train_data
  test = model_data$test_data
  ## Train Model
  ph = trainML_Models(trainLRs =train,
                      testLRs =test,
                      ytrain = ttData$y_train,
                      y_test = ttData$y_test,
                      testIDs = ttData$test_ids,
                      models = ensemble ) 
  ## Extract Performance and Predictions
  p = ph$performance[,c(1:3,15)];colnames(p)[4] = "model"
  preds= ph$predictionMatrix
  preds = data.frame(Seed = sd,
                     Fold = f,
                     Sparsity = max_sparsity,
                     num_features = number_parts,
                     num_samples = nsamps,
                     num_parts = ncol(model_data$train_data) ,
                     num_ratios = NA,
                     method = paste0("relAbundance_",methman),
                     comp_time = compTime[1],
                     preds)
  preds = left_join(preds,p)
  mdl = unique(preds$model)
  for(m in mdl){
    pr = preds %>% 
      filter(model==m)
    mroc = pROC::multiclass.roc(pr$Status,pr[,classes])
    c= combinat::combn2(classes)
    ## multiclass pairwise AUC
    if(length(classes)>2){
      ## Compute Pairwise AUC
      pwAUC = data.frame()
      for(jid in 1:nrow(c) ){
        pw = as.character(c[jid,])
        prn = pr
        prn = prn[prn$Status%in%pw,]
        mc = unique(prn$Status)
        controls = which(mc%in%c("nonIBD","UC"))
        cases = which(!mc%in%c("nonIBD","UC"))
        mc = c(mc[controls],mc[cases])
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
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = ncol(model_data$train_data) ,
                      num_ratios = NA,
                      method =  paste0("relAbundance_",methman),
                      comp_time = compTime[1],
                      p[p$model==m,],pairwise)
      perf.df = rbind(perf.df,pp)
      
    }else{
      ## Compute Optimal Threshold and Metrics
      optThreshold = pROC::roc(pr$Status,pr[,classes[1]],levels = classes)
      thrhold = pROC::coords(optThreshold,x="best",
                             ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                             transpose = FALSE)
      optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
      kk = which.max(optMCC)
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = ncol(model_data$train_data) ,
                      num_ratios = NA,
                      method = paste0("relAbundance_",methman),
                      comp_time = compTime[1],mcc = optMCC[kk],thrhold[kk,],
                      p[p$model==m,],AUC = as.numeric(pROC::auc(mroc)))
      perf.df = rbind(perf.df,pp)
    }
  }
  
}))

## ALR
message("Boruta - Process ALR")
suppressMessages(suppressWarnings({
  compTime = system.time({
    model_data = DiCoVarML::alrFeatureSelection(train_data = ttData$train_Data,
                                                y_train = ttData$y_train,num_borutaRuns = 200,
                                                test_data = ttData$test_data,
                                                featureSelectionMethod = selection_method,
                                                impute_factor = ttData$imputeFactor)
  })
  tr = model_data$train_data;colnames(tr) = paste0(colnames(tr),"___",colnames(ttData$train_Data)[1])
  ts = model_data$test_data;colnames(ts) = paste0(colnames(ts),"___",colnames(ttData$train_Data)[1])
  ## Train Model
  ph = trainML_Models(trainLRs = tr,
                      testLRs =ts,
                      ytrain = ttData$y_train,
                      y_test = ttData$y_test,testIDs = ttData$test_ids,models = ensemble ) 
  ## Extract Performance and Predictions
  cn = colnames(tr)
  n_ratios = length(cn)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)
  p = ph$performance[,c(1:3,15)];colnames(p)[4] = "model"
  
  ## Predictions
  preds= ph$predictionMatrix
  preds = data.frame(Seed = sd,
                     Fold = f,
                     Sparsity = max_sparsity,
                     num_features = number_parts,
                     num_samples = nsamps,
                     num_parts = length(uniqueParts) ,
                     num_ratios = n_ratios,
                     method = paste0("alr_",methman),
                     comp_time = compTime[1],
                     preds)
  preds = left_join(preds,p)
  mdl = unique(preds$model)
  for(m in mdl){
    pr = preds %>% 
      filter(model==m)
    mroc = pROC::multiclass.roc(pr$Status,pr[,classes])
    
    
    c= combinat::combn2(classes)
    ## multiclass pairwise AUC
    if(length(classes)>2){
      ## Compute Pairwise AUC
      pwAUC = data.frame()
      for(jid in 1:nrow(c) ){
        pw = as.character(c[jid,])
        prn = pr
        prn = prn[prn$Status%in%pw,]
        mc = unique(prn$Status)
        controls = which(mc%in%c("nonIBD","UC"))
        cases = which(!mc%in%c("nonIBD","UC"))
        mc = c(mc[controls],mc[cases])
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
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = paste0("alr_",methman),
                      comp_time = compTime[1],
                      p[p$model==m,],pairwise)
      perf.df = rbind(perf.df,pp)
      
    }else{
      ## Compute Optimal Threshold and Metrics
      optThreshold = pROC::roc(pr$Status,pr[,classes[1]],levels = classes)
      thrhold = pROC::coords(optThreshold,x="best",
                             ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                             transpose = FALSE)
      optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
      kk = which.max(optMCC)
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = paste0("alr_",methman),
                      comp_time = compTime[1],mcc = optMCC[kk],thrhold[kk,],
                      p[p$model==m,],AUC = as.numeric(pROC::auc(mroc)))
      perf.df = rbind(perf.df,pp)
    }
    
    
    
    
  }
}))

## PLR
message("Boruta - Process PLR")
suppressMessages(suppressWarnings({
  compTime = system.time({
    model_data = DiCoVarML::plrFeatureSelection(train_data = ttData$train_Data,
                                                y_train = ttData$y_train,num_borutaRuns = 200,
                                                test_data = ttData$test_data,
                                                featureSelectionMethod = selection_method,
                                                impute_factor = ttData$imputeFactor)
  })
  
  ## Train Model
  plr_featsBoruta = model_data
  ph = trainML_Models(trainLRs = model_data$train_data,
                      testLRs =model_data$test_data,
                      ytrain = ttData$y_train,
                      y_test = ttData$y_test,
                      testIDs = ttData$test_ids,
                      models = ensemble ) 
  ## Extract Performance and Predictions
  cn = colnames(model_data$train_data)
  n_ratios = length(cn)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)
  p = ph$performance[,c(1:3,15)];colnames(p)[4] = "model"
  ## Predictions
  preds= ph$predictionMatrix
  preds = data.frame(Seed = sd,
                     Fold = f,
                     Sparsity = max_sparsity,
                     num_features = number_parts,
                     num_samples = nsamps,
                     num_parts = length(uniqueParts) ,
                     num_ratios = n_ratios,
                     method = paste0("plr_",methman),
                     comp_time = compTime[1],
                     preds)
  preds = left_join(preds,p)
  preds_plr = preds
  mdl = unique(preds$model)
  for(m in mdl){
    pr = preds %>% 
      filter(model==m)
    mroc = pROC::multiclass.roc(pr$Status,pr[,classes])
    
    
    c= combinat::combn2(classes)
    ## multiclass pairwise AUC
    if(length(classes)>2){
      ## Compute Pairwise AUC
      pwAUC = data.frame()
      for(jid in 1:nrow(c) ){
        pw = as.character(c[jid,])
        prn = pr
        prn = prn[prn$Status%in%pw,]
        mc = unique(prn$Status)
        controls = which(mc%in%c("nonIBD","UC"))
        cases = which(!mc%in%c("nonIBD","UC"))
        mc = c(mc[controls],mc[cases])
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
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = paste0("plr_",methman),
                      comp_time = compTime[1],
                      p[p$model==m,],pairwise)
      perf.df = rbind(perf.df,pp)
      
    }else{
      ## Compute Optimal Threshold and Metrics
      optThreshold = pROC::roc(pr$Status,pr[,classes[1]],levels = classes)
      thrhold = pROC::coords(optThreshold,x="best",
                             ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                             transpose = FALSE)
      optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
      kk = which.max(optMCC)
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = paste0("plr_",methman),
                      comp_time = compTime[1],mcc = optMCC[kk],thrhold[kk,],
                      p[p$model==m,],AUC = as.numeric(pROC::auc(mroc)))
      perf.df = rbind(perf.df,pp)
    }
    
    
  }
}))



# DCV ---------------------------------------------------------------
dcv_feats = list()
train_control <- caret::trainControl(method="repeatedcv",
                                     repeats = 5,
                                     number=5,
                                     seeds = NULL,
                                     classProbs = TRUE,
                                     savePredictions = F,
                                     allowParallel = TRUE,
                                     summaryFunction = caret::multiClassSummary
)
## Compute DCV Filtered Features
compTime = system.time({
  feats = DiCoVarML::dcvRatioFilter(xtrain = ttData$train_Data,useKFN = T,
                                    th_percent = 0.25, 
                                    ytrain = ttData$y_train,
                                    xtest = ttData$test_data, 
                                    impute_factor = ttData$imputeFactor)
})

message("Fit MST-RFE subset")
## RFE
{
  compTime2 = system.time({
    pp = rfeSelection.ByMetric(train_ratio = feats$MST$train,
                               test_ratio = feats$MST$test,
                               ytrain =ttData$y_train,
                               ntrees = 750,
                               sets = 20,
                               impMeasure = "impurity",
                               kfold = 10,
                               minPercentFeatReturn = .05)
  })
  train_data2 = pp$reducedTrainRatios
  test_data2 = pp$reducedTestRatio
  dcv_feats[["mst_rfe"]] = list(train = train_data2,test =test_data2)
  message("number of features = ",ncol(train_data2))
  ## Train Model
  ph = trainML_Models(trainLRs = train_data2,
                      testLRs = test_data2,
                      ytrain = ttData$y_train,
                      y_test = ttData$y_test,
                      testIDs = ttData$test_ids,
                      models = ensemble ) 
  ## Extract Performance and Predictions
  cn = colnames(train_data2)
  n_ratios = length(cn)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)
  p = ph$performance[,c(1:3,15)];colnames(p)[4] = "model"
  
  ## Predictions
  preds= ph$predictionMatrix
  preds = data.frame(Seed = sd,
                     Fold = f,
                     Sparsity = max_sparsity,
                     num_features = number_parts,
                     num_samples = nsamps,
                     num_parts = length(uniqueParts) ,
                     num_ratios = n_ratios,
                     method = "dcv_mstRFE",
                     comp_time = compTime[1]+compTime2[1],
                     preds)
  preds = left_join(preds,p)
  mdl = unique(preds$model)
  for(m in mdl){
    pr = preds %>% 
      filter(model==m)
    mroc = pROC::multiclass.roc(pr$Status,pr[,classes])
    
    c= combinat::combn2(classes)
    ## multiclass pairwise AUC
    if(length(classes)>2){
      ## Compute Pairwise AUC
      pwAUC = data.frame()
      for(jid in 1:nrow(c) ){
        pw = as.character(c[jid,])
        prn = pr
        prn = prn[prn$Status%in%pw,]
        mc = unique(prn$Status)
        controls = which(mc%in%c("nonIBD","UC"))
        cases = which(!mc%in%c("nonIBD","UC"))
        mc = c(mc[controls],mc[cases])
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
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = "dcv_mstRFE",
                      comp_time = compTime[1]+compTime2[1],
                      p[p$model==m,],pairwise)
      
      perf.df = rbind(perf.df,pp)
      
    }else{
      ## Compute Optimal Threshold and Metrics
      optThreshold = pROC::roc(pr$Status,pr[,classes[1]],levels = classes)
      thrhold = pROC::coords(optThreshold,x="best",
                             ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                             transpose = FALSE)
      optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
      kk = which.max(optMCC)
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = "dcv_mstRFE",
                      comp_time = compTime[1]+compTime2[1],mcc = optMCC[kk],thrhold[kk,],
                      p[p$model==m,],
                      AUC = as.numeric(pROC::auc(mroc))
      )
      perf.df = rbind(perf.df,pp)
    }
    
  }
  preds_mst = preds
}

message("Fit MST-GLM subset")
## penalized regression
{
  compTime2 = system.time({
    glm.mdl1 = caret::train(x = as.matrix(feats$MST$train) ,
                            y =ttData$y_train,
                            metric = "ROC",
                            max.depth = 0,
                            method = "glmnet",
                            trControl = train_control
    )
  })
  
  imp  = caret::varImp(glm.mdl1)
  imp  = data.frame(feature = rownames(imp$importance),imp = imp$importance,total = rowSums(imp$importance))
  keep = imp[imp$total>0,]
  keep = keep$feature
  if(length(keep)>2){
    train_data2 = subset(feats$MST$train,select = c(keep))
    test_data2 = subset(feats$MST$test,select = c(keep))
  }else{
    train_data2 = feats$MST$train
    test_data2 =  feats$MST$test
  }
  dcv_feats[["mst_glm"]] = list(train = train_data2,test =test_data2)
  
  ## Train Model
  ph = trainML_Models(trainLRs = train_data2,
                      testLRs = test_data2,
                      ytrain = ttData$y_train,
                      y_test = ttData$y_test,
                      testIDs = ttData$test_ids,
                      models = ensemble ) 
  ## Extract Performance and Predictions
  cn = colnames(train_data2)
  n_ratios = length(cn)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)
  p = ph$performance[,c(1:3,15)];colnames(p)[4] = "model"
  
  
  ## Predictions
  preds= ph$predictionMatrix
  preds = data.frame(Seed = sd,
                     Fold = f,
                     Sparsity = max_sparsity,
                     num_features = number_parts,
                     num_samples = nsamps,
                     num_parts = length(uniqueParts) ,
                     num_ratios = n_ratios,
                     method = "dcv_mstGLM",
                     comp_time = compTime[1]+compTime2[1],
                     preds)
  preds = left_join(preds,p)
  mdl = unique(preds$model)
  for(m in mdl){
    pr = preds %>% 
      filter(model==m)
    mroc = pROC::multiclass.roc(pr$Status,pr[,classes])
    
    c= combinat::combn2(classes)
    ## multiclass pairwise AUC
    if(length(classes)>2){
      ## Compute Pairwise AUC
      pwAUC = data.frame()
      for(jid in 1:nrow(c) ){
        pw = as.character(c[jid,])
        prn = pr
        prn = prn[prn$Status%in%pw,]
        mc = unique(prn$Status)
        controls = which(mc%in%c("nonIBD","UC"))
        cases = which(!mc%in%c("nonIBD","UC"))
        mc = c(mc[controls],mc[cases])
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
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = "dcv_mstGLM",
                      comp_time = compTime[1]+compTime2[1],
                      p[p$model==m,],pairwise)
      
      perf.df = rbind(perf.df,pp)
      
    }else{
      ## Compute Optimal Threshold and Metrics
      optThreshold = pROC::roc(pr$Status,pr[,classes[1]],levels = classes)
      thrhold = pROC::coords(optThreshold,x="best",
                             ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                             transpose = FALSE)
      optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
      kk = which.max(optMCC)
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = "dcv_mstGLM",
                      comp_time = compTime[1]+compTime2[1],mcc = optMCC[kk],thrhold[kk,],
                      p[p$model==m,],
                      AUC = as.numeric(pROC::auc(mroc))
      )
      perf.df = rbind(perf.df,pp)
    }
    
  }
  preds_mstGLM = preds
}

message("Fit Dense-RFE subset")
## RFE
{
  compTime2 = system.time({
    pp = rfeSelection.ByMetric(train_ratio = feats$Dense$train,
                               test_ratio = feats$Dense$test,
                               ytrain =ttData$y_train,
                               ntrees = 750,
                               sets = 20,
                               impMeasure = "impurity",
                               kfold = 5,
                               minPercentFeatReturn = .05)
  })
  train_data2 = pp$reducedTrainRatios
  test_data2 = pp$reducedTestRatio
  dcv_feats[["dense_rfe"]] = list(train = train_data2,test =test_data2)
  message("number of features = ",ncol(train_data2))
  
  ## Train Model
  ph = trainML_Models(trainLRs = train_data2,
                      testLRs = test_data2,
                      ytrain = ttData$y_train,
                      y_test = ttData$y_test,
                      testIDs = ttData$test_ids,
                      models = ensemble ) 
  ## Extract Performance and Predictions
  cn = colnames(train_data2)
  n_ratios = length(cn)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)
  p = ph$performance[,c(1:3,15)];colnames(p)[4] = "model"
  ## Predictions
  preds= ph$predictionMatrix
  preds = data.frame(Seed = sd,
                     Fold = f,
                     Sparsity = max_sparsity,
                     num_features = number_parts,
                     num_samples = nsamps,
                     num_parts = length(uniqueParts) ,
                     num_ratios = n_ratios,
                     method = "dcv_denseRFE",
                     comp_time = compTime[1]+compTime2[1],
                     preds)
  preds = left_join(preds,p)
  mdl = unique(preds$model)
  for(m in mdl){
    pr = preds %>% 
      filter(model==m)
    mroc = pROC::multiclass.roc(pr$Status,pr[,classes])
    
    c= combinat::combn2(classes)
    ## multiclass pairwise AUC
    if(length(classes)>2){
      ## Compute Pairwise AUC
      pwAUC = data.frame()
      for(jid in 1:nrow(c) ){
        pw = as.character(c[jid,])
        prn = pr
        prn = prn[prn$Status%in%pw,]
        mc = unique(prn$Status)
        controls = which(mc%in%c("nonIBD","UC"))
        cases = which(!mc%in%c("nonIBD","UC"))
        mc = c(mc[controls],mc[cases])
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
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method = "dcv_denseRFE",
                      comp_time = compTime[1]+compTime2[1],
                      p[p$model==m,],pairwise)
      perf.df = rbind(perf.df,pp)
      
    }else{
      ## Compute Optimal Threshold and Metrics
      optThreshold = pROC::roc(pr$Status,pr[,classes[1]],levels = classes)
      thrhold = pROC::coords(optThreshold,x="best",
                             ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                             transpose = FALSE)
      optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
      kk = which.max(optMCC)
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method =  "dcv_denseRFE",
                      comp_time = compTime[1]+compTime2[1],mcc = optMCC[kk],thrhold[kk,],
                      p[p$model==m,],
                      AUC = as.numeric(pROC::auc(mroc))
      )
      perf.df = rbind(perf.df,pp)
    }
    
    
  }
  preds_dense = preds
}

message("Fit Dense-GLM subset")
## penalized regression
{
  compTime2 = system.time({
    glm.mdl1 = caret::train(x = as.matrix(feats$Dense$train) ,
                            y =ttData$y_train,
                            metric = "ROC",
                            max.depth = 0,
                            method = "glmnet",
                            trControl = train_control)
  })
  imp  = caret::varImp(glm.mdl1)
  imp  = data.frame(feature = rownames(imp$importance),imp = imp$importance,total = rowSums(imp$importance))
  keep = imp[imp$total>0,]
  keep = keep$feature
  if(length(keep)>2){
    train_data2 = subset(feats$Dense$train,select = c(keep))
    test_data2 = subset(feats$Dense$test,select = c(keep))
  }else{
    train_data2 = feats$Dense$train
    test_data2 =  feats$Dense$test
  }
  dcv_feats[["dense_glm"]] = list(train = train_data2,test =test_data2)
  
  ph = trainML_Models(trainLRs = train_data2,
                      testLRs = test_data2,
                      ytrain = ttData$y_train,
                      y_test = ttData$y_test,
                      testIDs = ttData$test_ids,
                      models = ensemble ) 
  ## Extract Performance and Predictions
  cn = colnames(train_data2)
  n_ratios = length(cn)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)
  p = ph$performance[,c(1:3,15)];colnames(p)[4] = "model"
  ## Predictions
  preds= ph$predictionMatrix
  preds = data.frame(Seed = sd,
                     Fold = f,
                     Sparsity = max_sparsity,
                     num_features = number_parts,
                     num_samples = nsamps,
                     num_parts = length(uniqueParts) ,
                     num_ratios = n_ratios,
                     method =  "dcv_denseGLM",
                     comp_time = compTime[1]+compTime2[1],
                     preds)
  preds = left_join(preds,p)
  mdl = unique(preds$model)
  for(m in mdl){
    pr = preds %>% 
      filter(model==m)
    mroc = pROC::multiclass.roc(pr$Status,pr[,classes])
    
    c= combinat::combn2(classes)
    ## multiclass pairwise AUC
    if(length(classes)>2){
      ## Compute Pairwise AUC
      pwAUC = data.frame()
      for(jid in 1:nrow(c) ){
        pw = as.character(c[jid,])
        prn = pr
        prn = prn[prn$Status%in%pw,]
        mc = unique(prn$Status)
        controls = which(mc%in%c("nonIBD","UC"))
        cases = which(!mc%in%c("nonIBD","UC"))
        mc = c(mc[controls],mc[cases])
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
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method =  "dcv_denseGLM",
                      comp_time = compTime[1]+compTime2[1],
                      p[p$model==m,],pairwise)
      
      perf.df = rbind(perf.df,pp)
      
    }else{
      ## Compute Optimal Threshold and Metrics
      optThreshold = pROC::roc(pr$Status,pr[,classes[1]],levels = classes)
      thrhold = pROC::coords(optThreshold,x="best",
                             ret=c("threshold","tn","tp","fp","fn","ppv","npv","sen","spe","youden"),
                             transpose = FALSE)
      optMCC = mltools::mcc(TP = thrhold$tp,FP = thrhold$fp,TN = thrhold$tn,FN = thrhold$fn)
      kk = which.max(optMCC)
      pp = data.frame(Seed = sd,
                      Fold = f,
                      Sparsity = max_sparsity,
                      num_features = number_parts,
                      num_samples = nsamps,
                      num_parts = length(uniqueParts) ,
                      num_ratios = n_ratios,
                      method =  "dcv_denseGLM",
                      comp_time = compTime[1]+compTime2[1],mcc = optMCC[kk],thrhold[kk,],
                      p[p$model==m,],
                      AUC = as.numeric(pROC::auc(mroc))
      )
      perf.df = rbind(perf.df,pp)
    }
    
    
  }
  preds_denseGLM = preds
}



# Write Output --------------------------------------------------------
out = list(predictions = list(mst = preds_mst,dense = preds_dense),
           features = list(DCV = feats),
           sample_info = list(train = ttData$train_ids,test = ttData$test_ids),
           labels = list(train = ttData$y_train,test = ttData$y_test))
fs = paste0("/nas/longleaf/home/andrew84/rProjects/dicovarFS/Fold",f,"_Seed",sd,"_Sparisty",max_sparsity*100,"_Dataset",
            dataset,"_selectionMethod",selection_method,"_featSets.Rda")
save(out,file = fs)
readr::write_csv(perf.df,fn)


