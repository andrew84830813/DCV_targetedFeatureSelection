


#' Targeted DCV Feature Selection
#'
#' @param trainx train data partition samples x parts/taxa/etc.
#' @param testx test data partition samples x parts/taxa/etc.
#' @param y_label train labels
#' @param y_test test labels
#' @param ensemble ensemble model definition
#' @param tarFeatures number of parts/taxa/etc. to retain (independent from number of ratios to retain)
#' @param imp_factor factor to multiplicative impute count zeros
#' @param ts.id test set id matrix. Can be found in output list from partitioning functions i.e.  kfoldDataPartition(),lodo_partition(),etc. Can be manually entered nrows = nrows(testx)
#' @param seed Random Seed control for reproducibility
#' @param max_sparsity max sparsity of parts/taxa/etc. For example max_sparsity=0.10 would mean to only retain parts/taxa/etc present in at least 10% of all samples.
#' @param select_randomFeatures should random parts/taxa/etc. be selected; useful for benchmarking
#' @param dcv a dcv matrix can optionally be provided
#' @param lrs.train An optional supplied logratio matrix for training data; must be supplied if dcv scores are pre computed
#' @param lrs.test An optional supplied logratio matrix for test data; must be supplied if dcv scores are pre computed
#' @param useRidgeWeights should ensemble model first be weighted by ridge regression coefficients?
#' @param scaledata should train data be scaled (test data scaled based on train)
#' @param use_rfe should a random forest recursive elimination model be trained 
#' @param minConnected should the ratio network be minimally connected or threshold with dcv_score>0
#'
#' @return
#' 
#' @export
#'
#' 
targeted_dcvSelection = function(trainx,
                                 testx,
                                 y_label,
                                 y_test,
                                 dcv=NULL,lrs.train=NULL,lrs.test = NULL,minConnected = F,
                                 ensemble  = c("ranger","pls","svmRadial","glmnet","rangerE"), tarFeatures = 5,
                                 imp_factor = 1e-7,
                                 select_randomFeatures = F,use_rfe = T,
                                 ts.id,seed = 08272008,max_sparsity = .9,useRidgeWeights=T,scaledata = T){
  
  result = data.frame()
  geo.mean = function(x){
    exp(mean(log(x)))
  }
  compTime = 0
  
  prob_matrices = list()
  
  
  # Train Final Model -------------------------------------------------------

  baseDims = ncol(trainx)
    
    if(!select_randomFeatures){
      
      ## Compute DCV
      message("Perform DCV Log Ratio Selection - (Step 1 of 4)")
      if(is.null(dcv)){
        compTime = system.time({
          
          ## impute
          trainx = data.frame(fastImputeZeroes(trainx,impFactor = imp_factor))
          testx = data.frame(fastImputeZeroes(testx,impFactor = imp_factor)) 
          
          ## compute log ratios
          lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = y_label,trainx))
          lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = y_test,testx))
          
          cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train, 
                                              includeInfoGain = T, nfolds = 1, numRepeats = 1, 
                                              rankOrder = F)
          
        }) 
        dcv = cc.dcv$lrs
      }
      
      ## Compute Node Strength
      dcv_strength = DiCoVarML::computeDCVStrength(list(dcv  =dcv))
      
      
      ## get subcomposition
      train_subcomp = subset(trainx,select = dcv_strength$Node[1:tarFeatures])
      test_subcomp = subset(testx,select = dcv_strength$Node[1:tarFeatures])
      
      
      ## get pair logratio from subcomp
      lrs.train = selEnergyPermR::calcLogRatio(df = data.frame(Status = y_label,train_subcomp))[,-1]
      lrs.test = selEnergyPermR::calcLogRatio(df = data.frame(Status = y_test,test_subcomp))[,-1]
      
      ## recompute dcvSCores
      if(minConnected){
        cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = data.frame(Status = y_label,lrs.train),
                                            includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                            rankOrder = T)
        dcv = cc.dcv$lrs
        w = min(which(dcv$nDistinct==tarFeatures))  
      }else{
        cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = data.frame(Status = y_label,lrs.train),
                                            includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                            rankOrder = F)
        dcv = cc.dcv$lrs
        w = max(which(dcv$rowmean>0))  
      }
      
      
      
      glm.train = subset(lrs.train,select = dcv$Ratio[1:w])
      glm.test = subset(lrs.test,select = dcv$Ratio[1:w])
      dcv = dcv[1:w,]
      
      
      
    }else{
      cn = colnames(trainx)
      compTime = 0
      rand_feats = sample(cn,replace = F)[1:tarFeatures]
      train_subcomp = subset(trainx,select = rand_feats)
      test_subcomp = subset(testx,select = rand_feats)
      ## get pair logratio from subcomp
      glm.train = selEnergyPermR::calcLogRatio(df = data.frame(Status = y_label,train_subcomp))[,-1]
      glm.test = selEnergyPermR::calcLogRatio(df = data.frame(Status = y_test,test_subcomp))[,-1]
    }
    
    
    
    
    ## Compute Number of Part
    cn = colnames(glm.train)
    uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
    n_parts  = dplyr::n_distinct(uniqueParts)
    
    
    ## scale data
    if(scaledata ){
      pp = caret::preProcess(glm.train,method = "scale")
      glm.train <- predict(pp, glm.train)
      glm.test     <- predict(pp, glm.test)
    }
    
    
    ## Train Ridge Model
    message("Train Ridge Regression Model - (Step 2 of 4)")
    
    type_family = dplyr::if_else(length(classes)>2,"multinomial","binomial")
    compTime2 = system.time({
      cv.clrlasso <- glmnet::cv.glmnet(as.matrix(glm.train),y_label, standardize=F, alpha=0,family=type_family)
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
    p = stats::predict(cv.clrlasso, newx = as.matrix(glm.test), s = "lambda.min",type = "response")
    model_ = list(mdl = cv.clrlasso,data = list(train = glm.train,test = glm.test))
    if(type_family=="binomial"){
      mroc = pROC::roc(y_test,p)
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
      pmat = data.frame(p,1-p);colnames(pmat) = classes

    }else{
      ## multiclass
      mroc = pROC::multiclass.roc(y_test,p[,,1])
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
      pmat = data.frame(p[,,1])#;colnames(pmat) = classes
      
    }
    ## Save Performance
    perf = data.frame(Approach = "DCV-ridgeRegression",
                      AUC = as.numeric(mroc.dcvlasso),
                      number_parts = n_parts,
                      number_ratios = ncol(glm.train),
                      comp_time = compTime[3]+compTime2[3],
                      base_dims = baseDims)
    result = rbind(result,perf)
    prob_matrices[["ridgeRegression"]] = data.frame(Status = y_test,pmat)
    
    
    
    
    message("Train Ensemble Model - (Step 3 of 4)")
    ## Ridge Features with Ensemble  
    # use ridge weights
    if(type_family=="binomial"){
      train_data2 = subset(glm.train,select = unique(names(features)))
      test_data2 = subset(glm.test,select = unique(names(features)))
      
      
      if(useRidgeWeights){
        train_data2 = sweep(train_data2,MARGIN = 2,STATS = c,FUN = "*")
        test_data2 = sweep(test_data2,MARGIN = 2,STATS = c,FUN = "*")   
      }
       
    }else{
      train_data2 = subset(glm.train,select = feat.df$Ratio)
      test_data2 = subset(glm.test,select = feat.df$Ratio)
      
      if(useRidgeWeights){
        train_data2 = sweep(train_data2,MARGIN = 2,STATS = feat.df$coef,FUN = "*")
        test_data2 = sweep(test_data2,MARGIN = 2,STATS = feat.df$coef,FUN = "*")  
      }
      
    }
    
    weight.train = train_data2
    weight.test = test_data2
    
    ## Train Model
    ph = trainML_Models(trainLRs = train_data2,
                        testLRs = test_data2,
                        ytrain = y_label,
                        y_test = y_test,
                        testIDs = ts.id,
                        models = ensemble ) 
    pmat = ph$predictionMatrix
    pmat = pmat %>% 
      dplyr::group_by(ID,Status) %>% 
      dplyr::select(-model) %>% 
      dplyr::summarise_all(.funs = mean)
    pmat = data.frame(pmat)
    mroc = pROC::multiclass.roc(pmat$Status,pmat[,classes]);mroc
    ## Compute Number of Part
    cn = colnames(train_data2)
    n_ratios = length(cn)
    uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
    n_parts  = dplyr::n_distinct(uniqueParts)
    ## Save Performance
    perf = data.frame(Approach = "DCV-ridgeEnsemble",AUC = as.numeric(pROC::auc(mroc)),
                      number_parts = n_parts,number_ratios = ncol(train_data2),
                      comp_time = compTime[3]+compTime2[3],
                      base_dims = baseDims)
    result = rbind(result,perf)
    prob_matrices[["ridgeEnsemble"]] = pmat
    
    
    message("Train RFE Model - (Step 4 of 4)")
    rfe_Features = NULL
    if(use_rfe){
      suppressMessages(suppressWarnings({
        
        type = "Dense"
        ## Feature Selection: rf-RFE 
        tc = .99
        c1.cor = cor(glm.train,method = "spearman")
        c.fc = data.frame(Ratio = caret::findCorrelation(c1.cor,cutoff = tc,names = T))
        keep =!colnames(glm.train)%in%c.fc$Ratio
        glm.train = glm.train[,keep]
        glm.test = glm.test[,keep]
        
        compTime2 = system.time({
          pp = rfeSelection.ByMetric(train_ratio = glm.train,
                                     test_ratio = glm.test,
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
                            y_test = y_test,
                            testIDs = ts.id,
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
        
        ## Save Performance
        perf = data.frame(Approach = "DCV-rfRFE",AUC = as.numeric(pROC::auc(mroc)),
                          number_parts = n_parts,number_ratios = ncol(train_data2),
                          comp_time = compTime[3]+compTime2[3],
                          base_dims = baseDims)
        result = rbind(result,perf)
        prob_matrices[["rfRFE"]] = pmat
      }))
      rfe_Features = list(train = train_data2,test =test_data2)
    }
    
    
    
    
    
    return(list(Performance = result,all_model_preds = ph$predictionMatrix,
                glm_model = model_,
                weighted_features = list(train = weight.train,test = weight.test),
                part_matrices = list(train = train_subcomp,test = test_subcomp),
                ridge_coefficients = c,probMatrices  = prob_matrices,
                final_dcv = dcv, 
                ensemble_pmat = pmat,
                rfe_features = rfe_Features,
                ridge_pmat = p ))
    

  
  
}



