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
source(file = 'CODA_Functions/functions_coda_penalized_regression.R')



# Setup Cluster ---------------------------------------------------------------

## Detect Cores
# clus <- parallel::makeCluster(10)
# doParallel::registerDoParallel(clus)




# Read External Args ---------------------------------------------------------------

args = c(1,2,0)
args = commandArgs(trailingOnly = TRUE)
seed_ = as.numeric(args[1]) # random seed selection
shift_parm = as.numeric(args[2])
permute_labels = as.logical(as.numeric(args[3])) #should be 0(False) or 1(True)


## set random seed
set.seed(seed_)


## Scenario-1 = Sim from 16S Data with mean shift
g1 = 100
g2 = 100
sparsePercent = .9
benchmark = data.frame()




mdwgs = readRDS("Output/16sModel.RDS");
set.seed(seed_);
dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts));
dat = sample_n(dat,size = g1+g2,replace = F);
labels = sample(c(rep("S1",g1),rep("S2",g2)));
dat = data.frame(Status = labels,dat);
f_name = paste0("16S_meanShift",shift_parm,"_permute",permute_labels,"_seed",seed_);
#process shift
procData = processCompData(dat,minPrevalence = sparsePercent);
dat = procData$processedData;
impFact = procData$impFactor;
y = dat[,-1];
bool = colSums(y)==0;
y = y[,!bool];
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact));


sets = seq(1,1.3,length.out = 6)
df = simFromExpData.largeMeanShft(raMatrix = dat[,-1],
                                  n1 = g1,n2 = g2,
                                  featureShiftPercent =  sets[shift_parm],
                                  perFixedFeatures = .9)



for(sd in 1:5){
  set.seed(sd)
  
  ## Partition Data
  k_fold = 2
  overll_folds = caret::createFolds(df[,1],k =k_fold,list = F)
  allData = lodo_partition(df,dataset_labels = overll_folds,sd)
  
  message("\n\n``````` Start Seed:  ",sd,"````````````\n\n")
  
  #  Perform 2-fold cross validation -------------------------------
  
  for(f in 1:k_fold){
    
    
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
    
    # clr lasso ---------------------------------------------------------------
    
    suppressMessages(suppressWarnings({
      xt =train.data
      yt = ttData$y_train
      
      z <- log(xt+1)
      ztrain <- apply(z,2,function (x) x-rowMeans(z))
      
      z <- log(test.data+1)
      ztest <- apply(z,2,function (x) x-rowMeans(z))
      
      ##scale data
      pp = caret::preProcess(ztrain,method = "scale")
      ztrain <- predict(pp, ztrain)
      ztest     <- predict(pp, ztest)
      
      
      type_family = if_else(length(classes)>2,"multinomial","binomial")
      compTime = system.time({
        cv.clrlasso <- glmnet::cv.glmnet(ztrain, yt, standardize=F, alpha=1,family=type_family,parallel = T)
      })
      
      if(type_family=="binomial"){
        features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
        features = features[-1,]
        features = features[abs(features)>0]
        n_parts = length(features)
        
      }else{
        features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
        keep  = c()
        for(o in 1:length(features)){
          ph = as.matrix(features[[o]])
          feat = ph[-1,]
          feat = feat[abs(feat)>0]
          keep = c(keep,feat)
        }
        
        features = unique(keep)
        n_parts=length(features)
        
      }
      ## make predictions
      p = predict(cv.clrlasso, newx = ztest, s = "lambda.min",type = "response")
      if(type_family=="binomial"){
        mroc = pROC::roc(ttData$y_test,p)
        mroc.clrlasso = pROC::auc(mroc);mroc.clrlasso
      }else{
        ## multiclass
        mroc = pROC::multiclass.roc(ttData$y_test,p[,,1])
        mroc.clrlasso = pROC::auc(mroc);mroc.clrlasso
        
      }
      #Write Results
      perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                        Dataset = f_name,Seed = sd,Fold = f,Approach = "CLR-LASSO",
                        AUC = as.numeric(mroc.clrlasso),
                        number_parts = n_parts,number_ratios = n_parts ,comp_time = compTime[3],
                        base_dims = ncol(train.data))
      benchmark = rbind(benchmark,perf)
    }))
    
    
    # coda lasso --------------------------------------------------------------
    
    suppressMessages(suppressWarnings({
      xt =train.data
      yt = ttData$y_train
      ytr = as.numeric(yt)-1
      xt =train.data+1
      xtest = test.data+1
      
      ztransform = function(xt,p_c = NULL){
        ztrain <- log(xt)
        # z=matrix of covariates: add a first column of 1's for beta0
        ztrain <- cbind(rep(1,nrow(ztrain)),ztrain)
        ztrain <- as.matrix(ztrain)
        # c=linear constraint sum(betas)=0 (except beta0)
        c <- c(0,rep(1,ncol(xt)))
        c <- c/sqrt(normSqr(c))
        c <- as.matrix(c)
        
        if(is.null(p_c)){
          p_c <- c%*%t(c)
        }
        
        # p_cginv <- ginv(p_c);  #use this only if necessary
        p_c_cginv <- p_c; #p_c%*%p_cginv; #this product in our case is = p_c
        # z transformation (centering) for improvement of optimization
        # this transformation does not affect the estimation since the linear predictor is the same
        ztrain <- (ztrain-(ztrain%*%p_c))
        colnames(ztrain) <- c("beta0",colnames(xt))
        return(list(dat = ztrain,proj = p_c))
      }
      
      ztrain = ztransform(xt)
      ztest = ztransform(xtest,p_c = ztrain$proj)$dat
      ztrain = ztrain$dat
      
      compTime = system.time({
        devexp = foreach::foreach(l = cv.clrlasso$lambda,.combine = rbind )%dopar%{
          HFHS.results_codalasso <- coda_logistic_lasso(ytr,(xt),lambda = l)
          de = HFHS.results_codalasso$`proportion of explained deviance`
          ph = data.frame(lambda = l,deviance_explained = de,num_variables = HFHS.results_codalasso$`number of selected variables`)
          ph$score = ph$deviance_explained *ph$lambda
          ph
        }
        devexp = devexp %>% 
          filter(num_variables>0)
        mx_dev = max(devexp$deviance_explained)
        devexp1 = devexp %>% 
          filter(deviance_explained == mx_dev) %>% 
          arrange(desc(lambda))
        
        HFHS.results_codalasso <- coda_logistic_lasso(ytr,(xt),lambda=devexp$lambda[1])
        n_parts = HFHS.results_codalasso$`number of selected variables`
        bet = HFHS.results_codalasso$betas
        train.balances = ztrain %*%bet
        test.balances = ztest %*%bet 
      })
      
      ## get train data
      df.train =  data.frame(bal = train.balances)
      df.test =  data.frame(bal = test.balances)
      
      
      tbl = data.frame(Status = ytr,df.train)
      tbl.test = data.frame(Status = ttData$y_test,df.test)
      gm = glm(formula = Status~.,data = tbl,family = binomial)
      p = predict.glm(gm, newdata = df.test, type = "response")
      mroc.codalasso = pROC::auc(ttData$y_test,p);mroc.codalasso
      #Write Results
      perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                        Dataset = f_name,Seed = sd,Fold = f,Approach = "Coda-LASSO",AUC = as.numeric(mroc.codalasso),
                        number_parts = n_parts,number_ratios = 1 ,comp_time = compTime[3],
                        base_dims = ncol(train.data)
      )
      benchmark = rbind(benchmark,perf)
    }))  
    
    
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
    
    
    
  }
}


benchmark$permuteLabel = permute_labels
benchmark$shift_parm = shift_parm
write_csv(x = benchmark,file = paste0("Results/",f_name,".csv"))
