
relAbunApproach = function(method){
  
  suppressMessages(suppressWarnings({
    model_data = DiCoVarML::relAbundanceFeatureSelection(train_data = ttData$train_Data,
                                                         y_train = ttData$y_train,
                                                       test_data = ttData$test_data,
                                                       featureSelectionMethod = 1,
                                                       impute_factor = ttData$imputeFactor )
  }))
  
  train_data = model_data$train_data
  test_data = model_data$test_data
  model_performance = DiCoVarML::processModel(train_x = train_data,train_y = ttData$y_train,
                                              test_x = test_data,
                                              test_y = ttData$y_test,
                                              test_ids = ttData$test_ids,
                                              ensembleModels = ensemble,
                                              train_stackedModel = trainStacked,
                                              avgModelPrediction = T)
  data.frame(Seed = sd,
             Fold = f,
             num_parts = ncol(train_data) ,
             num_ratios = NA,
             keyParms = NA,
             Approach = "rel_abund",
             Method = method,
             model_performance$performance)
}


plrApproach = function(method){
  
  suppressMessages(suppressWarnings({
    model_data = DiCoVarML::plrFeatureSelection(train_data = ttData$train_Data,y_train = ttData$y_train,
                                                test_data = ttData$test_data,
                                                featureSelectionMethod = method,
                                                impute_factor = ttData$imputeFactor)
  }))
  
  train_data = model_data$train_data
  test_data = model_data$test_data
  model_performance = DiCoVarML::processModel(train_x = train_data,train_y = ttData$y_train,
                                              test_x = test_data,test_y = ttData$y_test,
                                              test_ids = ttData$test_ids,
                                              train_stackedModel = trainStacked,
                                              ensembleModels = ensemble,
                                              avgModelPrediction = T)
  cn = colnames(train_data)
  n_ratios = length(cn)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)
  
  data.frame(Seed = sd,
             Fold = f,
             num_parts = n_parts ,
             num_ratios = n_ratios,
             keyParms = NA,
             Approach = "plr",
             Method = method,
             model_performance$performance)
  
}


alrApproach = function(method){
  
  suppressMessages(suppressWarnings({
    model_data = DiCoVarML::alrFeatureSelection(train_data = ttData$train_Data,y_train = ttData$y_train,
                                              test_data = ttData$test_data,
                                              featureSelectionMethod = method,
                                              impute_factor = ttData$imputeFactor)
  }))
  
  train_data = model_data$train_data
  test_data = model_data$test_data
  
  model_performance = DiCoVarML::processModel(train_x = train_data,train_y = as.character(ttData$y_train),
                                              test_x = test_data,test_y = ttData$y_test,
                                              test_ids = ttData$test_ids,
                                              train_stackedModel = trainStacked,
                                              ensembleModels = ensemble,
                                              avgModelPrediction = T)
  data.frame(Seed = sd,Fold = f,
               num_parts = ncol(train_data)+1 ,
               num_ratios = ncol(train_data),
               keyParms = NA,
               Approach = "alr",
               Method = method,
               model_performance$performance)
}








fitModel = function(sd,f,approach,method,ensemble,
                      train_data,test_data,y_train,
                      testIDs,num_parts,
                      y_test,train_avgEnsemble = T,
                      trainStacked = F){

  suppressMessages(suppressWarnings({
    ## Train Model
    model_performance = DiCoVarML::processModel(train_x = train_data,train_y = y_train,
                                                test_x = test_data,
                                                test_y = y_test,
                                                ensembleModels = ensemble,
                                                test_ids = testIDs,
                                                train_stackedModel = trainStacked,
                                                avgModelPrediction = train_avgEnsemble)
    ## Compute Number of Ratios and Parts
    cn = colnames(train_data)
    if(stringr::str_detect(cn[1],"___")){
      n_ratios = length(cn)
      uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
      n_parts  = dplyr::n_distinct(uniqueParts)
    }else{
      n_ratios = NA
      n_parts = ncol(train_data)
    }
  }))
  mp = model_performance$performance
  mp$AUC = as.numeric(mp$AUC)
  ## Store Performance
    list(performance = data.frame(Seed = sd,
                  Fold = f,
                  num_features = num_parts,
                  num_samples = nrow(train_data),
                  num_parts = n_parts ,
                  num_ratios = n_ratios,
                  Approach = approach,
                  Method = method,
                  mp),
              modelInfo = model_performance
            )
}






dcvApproach1 = function(target_features){
  
  ph.df = data.frame()
  for(key_parm in target_features ){
    ## Perform DCV Feature Selection (Discover in Train Data Only)
    model_data = DiCoVarML::dcvFeatureSelection(train_data = ttData$train_Data,
                                                y_train = ttData$y_train,
                                                test_data = ttData$test_data,
                                                dcv_matrix = dcv_mat$dcv,
                                                impute_factor = ttData$imputeFactor,
                                                rf_importance = "impurity_corrected",
                                                dcv_method = 1,num_sets = 10,
                                                tarFeats =key_parm,minFeats = key_parm)
    
    ## Extract Key Features
    train_data = model_data$train_data
    test_data = model_data$test_data
    
    ## Train Model
    model_performance = DiCoVarML::processModel(train_x = train_data,train_y = ttData$y_train,
                                                test_x = test_data,test_y = ttData$y_test,
                                                ensembleModels = ensemble,
                                                test_ids = ttData$test_ids,
                                                train_stackedModel = trainStacked,
                                                avgModelPrediction = T)
    ## Compute Number of Ratios and Parts
    cn = colnames(train_data)
    n_ratios = length(cn)
    uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
    n_parts  = dplyr::n_distinct(uniqueParts)
    
    ## Store Performance
    ph = data.frame(Seed = sd,
               Fold = f,
               num_parts = n_parts ,
               num_ratios = n_ratios,
               keyParms = key_parm,
               Approach = "dcv",
               Method = 1,
               model_performance$performance)
    
    ph.df = rbind(ph.df,ph)
  }
  
  ph.df
  
}


dcvApproach2 = function(upperbound_percentage){
  

  ## Perform DCV Feature Selection (Discover in Train Data Only)
  ph.df = data.frame()
  for(key_parm in upperbound_percentage ){
    model_data =  DiCoVarML::dcvFeatureSelection(train_data = ttData$train_Data,
                                                 y_train = ttData$y_train,
                                                 test_data = ttData$test_data, ## only required for feature return
                                                 dcv_matrix = dcv_mat$dcv,
                                                 impute_factor = ttData$imputeFactor,
                                                 rf_importance = "impurity",
                                                 dcv_method = 2,
                                                 tarFeats = round(key_parm*number_parts),
                                                 upperBound_percent = key_parm)
    
    ## Extract Key Features
    train_data = model_data$train_data
    test_data = model_data$test_data
    
    ## Train Model
    model_performance = DiCoVarML::processModel(train_x = train_data,train_y = ttData$y_train,
                                                test_x = test_data,test_y = ttData$y_test,
                                                ensembleModels = ensemble,
                                                test_ids = ttData$test_ids,
                                                train_stackedModel = trainStacked,
                                                avgModelPrediction = T)
    ## Compute Number of Ratios and Parts
    cn = colnames(train_data)
    n_ratios = length(cn)
    uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
    n_parts  = dplyr::n_distinct(uniqueParts)
    mp  =model_performance$performance
    
    
    ## Store Performance
    ph = data.frame(Seed = sd,
                    Fold = f,
                    num_parts = n_parts ,
                    num_ratios = n_ratios,
                    keyParms = key_parm,
                    Approach = "dcv",
                    Method = 2,
                    mp
                    )
    
    ph.df = rbind(ph.df,ph)
  }

  ph.df
}


dcvApproach3 = function(upperbound_percentage){
  ## Perform DCV Feature Selection (Discover in Train Data Only)
  ph.df = data.frame()
  for(key_parm in upperbound_percentage ){
    model_data =  DiCoVarML::dcvFeatureSelection(train_data = ttData$train_Data,
                                                 y_train = ttData$y_train,
                                                 test_data = ttData$test_data, ## only required for feature return
                                                 dcv_matrix = dcv_mat$dcv,
                                                 impute_factor = ttData$imputeFactor,
                                                 rf_importance = "impurity",
                                                 dcv_method = 3,
                                                 upperBound_percent = key_parm)
    
    ## Extract Key Features
    train_data = model_data$train_data
    test_data = model_data$test_data
    
    ## Train Model
    model_performance = DiCoVarML::processModel(train_x = train_data,train_y = ttData$y_train,
                                                test_x = test_data,test_y = ttData$y_test,
                                                ensembleModels = ensemble,
                                                test_ids = ttData$test_ids,
                                                train_stackedModel = trainStacked,
                                                avgModelPrediction = T)
    ## Compute Number of Ratios and Parts
    cn = colnames(train_data)
    n_ratios = length(cn)
    uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
    n_parts  = dplyr::n_distinct(uniqueParts)
    
    ## Store Performance
    ph = data.frame(Seed = sd,
                    Fold = f,
                    num_parts = n_parts ,
                    num_ratios = n_ratios,
                    keyParms = key_parm,
                    Approach = "dcv",
                    Method = 3,
                    model_performance$performance)
    
    ph.df = rbind(ph.df,ph)
  }
  
  ph.df

}
