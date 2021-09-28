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

# # Detect Cores
# clus <- parallel::makeCluster(10)
# doParallel::registerDoParallel(clus)




# Read External Args ---------------------------------------------------------------

args = c(2,4,0)
args = commandArgs(trailingOnly = TRUE)
sd = as.numeric(args[1]) # random seed selection
dataset = as.numeric(args[2])
permute_labels = as.logical(as.numeric(args[3])) #should be 0(False) or 1(True)





# Load Data ---------------------------------------------------------------

switch(dataset,
       
       {
         load("C:/Users/andrew84/Documents/MicrobiomeProject/Data/HIV.rda");
         df = data.frame(Status = HIV$HIV_Status,HIV[,1:60]);
         f_name = "selbal_HIV_16s";
       },
       
       {
         load("C:/Users/andrew84/Documents/MicrobiomeProject/Data/Crohn.rda")
         df = data.frame(Status = Crohn$y,Crohn[,1:48]);
         df$Status = dplyr::if_else(df$Status=="CD","pos","neg");
         #preprocessing
         #remove g__ and f__ in taxa names
         taxa = data.frame(Taxa = colnames(df[,-1]));
         taxa$taxa = stringr::str_replace(taxa$Taxa,pattern = "\\.",replace = "");
         taxa = tidyr::separate(taxa,col = 1,into = c("des","taxa"),sep = "__");
         colnames(df) = c("Status",taxa$taxa);
         rs= rowSums(df[,-1])
         bool = rs>1e3;
         df = df[bool,];
         f_name = "selbal_Crohns_16s";
       },
       
       {
         ## nfald
         df = read.csv(file = "/Output/16S-11635_NAFLD.csv") ## good;
         rs = rowSums(df[,-1]);
         bool = rs>5e3;
         df = df[bool,];
         f_name = "qitta_NAFLD_16s";
         
       },
       
       {
         load("Output/microbiomeHD_cdiSchubert.Rda"); # a go
         df = expResults$taxa.df;
         df = data.frame(sample_id = rownames(df),df);
         df$Status = str_replace_all(df$Status,pattern = "\\/","_");
         df = df[!is.na(df$Status),];
         keep = unique(df$Status);
         keep = c("H","CDI");
         df = df[df$Status%in%keep,];
         df = df[,-1];
         f_name = "mbiomeHD_cdiSchubert_16s";
         
       },
       
       {
         load("Output/curatedMetaGenome_RubelMA_2020.Rda"); ## good
         keep = c("control","STH");
         f_name = "cmg_RubelMA-2020_STH"
       },
       
       {
         load("Output/curatedMetaGenome_ZhuF_2020.Rda") ## good
         keep = c("control","schizofrenia");
         f_name = "cmg_ZhuF-2020_schizo"
         
       },
       
       {
         load("Output/curatedMetaGenome_QinN_2014.Rda")
         keep = c("control","cirrhosis");
         f_name = "cmg_QinN-2014_cirr"
       },
       
       {
         load("Output/curatedMetaGenome_NielsenHB_2014.Rda");
         keep = c("control","IBD");
         f_name = "cmg_NielsenHB-2014_ibd"
         
       },
       
       {
         load("Output/curatedMetaGenome_FengQ_2015.Rda");
         keep = c("control","CRC");
         f_name = "cmg_FengQ-2015_crc"
         
       },
       
       {
         load("Output/curatedMetaGenome_ThomasAM_2019_c.Rda");
         keep = c("control","CRC");
         f_name = "cmg_ThomasAM_2019_crc"
         
       },
       
       {
         load("Output/curatedMetaGenome_WirbelJ_2018.Rda");
         keep = c("control","CRC");
         f_name = "cmg_WirbelJ-2018_crc"
       },
       
       {
         load("Output/curatedMetaGenome_YachidaS_2019.Rda");
         f_name = "cmg_YachidaS-2019_crc";
         keep = c("control","CRC")       
       },
       
       {
         load("Output/curatedMetaGenome_ZellerG_2014.Rda");
         f_name = "cmg_ZellerG_2014_crc";
         keep = c("control","CRC")
       },
       
       {
         load("Output/curatedMetaGenome_VogtmannE_2016.Rda");
         f_name = "cmg_ZVogtmannE_2016_crc";
         keep = c("control","CRC")
       },
       
       {
         load("Output/curatedMetaGenome_YuJ_2015.Rda");
         f_name = "cmg_YuJ_2015_crc";
         keep = c("control","CRC")
       }
       
)






if(str_detect(f_name,"cmg")){
  ## Process Data
  df = expResults$taxa.df
  df = data.frame(sample_id = rownames(df),df)
  df$Status = str_replace_all(df$Status,pattern = "\\/","_")
  rs = rowSums(df[,-2:-1])
  bool = rs>10e6
  df = df[bool,]
  df = df[df$Status %in% keep,]
  md = expResults$metadata
  ## by taxa level
  cname.table = expResults$taxa_info [,-1]
  taxa_levels =colnames(cname.table)
  cname.table = data.frame(ID = cname.table$species,cname.table)
  rdphf = tidyr::gather(df,"ID","Count",3:ncol(df))
  cname.table$ID = stringr::str_replace_all(cname.table$ID,pattern = " ",replacement = "\\.")
  rdphf$ID = stringr::str_replace_all(rdphf$ID,pattern = " ",replacement = "\\.")
  rdphf = left_join(rdphf,cname.table)
  level.list = list()
  for(taxa_level in c("species")){
    rdf = subset(rdphf,select = c("sample_id","Status",taxa_level,"Count"))
    rdf = rdf %>%
      dplyr::group_by(sample_id,get(taxa_level),Status) %>%
      dplyr::summarise(count = sum(Count) ) %>%
      rename(tlevel = "get(taxa_level)" ) %>%
      spread("tlevel","count",fill = 0)
    rdf = data.frame(rdf[,-1],row.names = rdf$sample_id )
    # rdf = rdf[rdf$Status %in% c("control","CRC"),]
    colnames(rdf)[-1] = paste0(str_sub(taxa_level,1,1),"__",colnames(rdf)[-1])
    level.list[[taxa_level]] = rdf
  }
  df = level.list$species
  df = df[!is.na(df$Status),]
}


f_name = paste0(f_name,"_seed",sd)
mxSparsePercent = .9
benchmark = data.frame()
seed_ = sd

#for(sd in 1:5){
set.seed(sd)

## Partition Data
k_fold = 2
overll_folds = caret::createFolds(df[,1],k =k_fold,list = F)
allData = lodo_partition(df,dataset_labels = overll_folds,sd)

## DCV Parms
scale_data = T
performRFE = F
useRidgeWeight = F
min_connected = F
ensemble = c("ranger","xgbTree","xgbLinear")
ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")



message("\n\n``````` Start Seed:  ",sd,"````````````\n\n")

#  Perform 2-fold cross validation -------------------------------

for(f in 1:k_fold){
  
  
  ## Extract Test/Train Split
  ttData = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                            fold = f,
                                            permLabels = permute_labels,
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
  
  
  
  # # selbal ------------------------------------------------------------------
  
  xt =train.data
  yt = ttData$y_train
  
  
  
  compTime = system.time({
    bv = selbal::selbal.aux(x = xt, y = yt,zero.rep = "one")
    train.bv = selbal::bal.value(bv,x = log(xt+1) )
  })

  
  df.train =  data.frame(bal = train.bv)
  test.bv = selbal::bal.value(bv,x = log(ttData$test_data +1))
  df.test =  data.frame(bal = test.bv)
  
  
  tbl = data.frame(Status = yt,df.train)
  tbl.test = data.frame(Status = ttData$y_test,df.test)
  gm = glm(formula = Status~.,data = tbl,family = binomial)
  p = predict.glm(gm, newdata = df.test, type = "response")
  mroc.selbal = pROC::auc(ttData$y_test,p);mroc.selbal
  
  
  
  #Write Results
  perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                    Dataset = f_name,Seed = sd,Fold = f,Approach = "SELBAL",AUC = as.numeric(mroc.selbal),
                    number_parts = length(bv$Taxa),number_ratios = 1 ,comp_time = compTime[3],
                    base_dims = ncol(train.data)
  )
  
  benchmark = rbind(benchmark,perf)
  
 
  
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
  
  
  #ensemble = c("ranger","xgbTree","xgbLinear")
  max_sparsity = .9
  
  
  ## Tune target features
  compTime1 = system.time({
    
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
          
          tar_dcvInner = targeted_dcvSelection(trainx = trainx,minConnected = min_connected,
                                               useRidgeWeights = useRidgeWeight,use_rfe = performRFE,scaledata = scale_data,
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
  
  })
    
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
    
    
    ## Apply targted feature selection method
    compTime2 = system.time({
      cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train, 
                                          includeInfoGain = T, nfolds = 1, numRepeats = 1, 
                                          rankOrder = F)
      tar_Features = inner_perf2$tar_Features[which.max(inner_perf2$AUC)]
      tar_dcv = targeted_dcvSelection(trainx = trainx,minConnected = min_connected,
                                      useRidgeWeights = useRidgeWeight,use_rfe = performRFE,scaledata = scale_data,
                                      testx = testx,
                                      dcv = cc.dcv$lrs,lrs.train = lrs.train,lrs.test = lrs.test,
                                      y_label = ttData$y_train,
                                      seed = sd,
                                      #ensemble = ensemble,
                                      y_test = ttData$y_test,
                                      tarFeatures = tar_Features,
                                      ts.id = ttData$test_ids, 
                                      max_sparsity = max_sparsity
      )
  })
  
  perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                    Dataset = f_name,Seed = sd,Fold = f,Approach = tar_dcv$Performance$Approach,AUC = as.numeric(tar_dcv$Performance$AUC),
                    number_parts = tar_dcv$Performance$number_parts,number_ratios = tar_dcv$Performance$number_ratios ,
                    comp_time = compTime1[3]+compTime2[3],
                    base_dims = ncol(train.data)
  )
  benchmark = rbind(benchmark,perf)
  
  
  
}

#}



write_csv(x = benchmark,file = paste0("Results/Exp_Benchamark/",f_name,".csv"))


# res = benchmark %>% 
#   group_by(Approach,Dataset) %>% 
#   summarise_all(mean) %>% 
#   arrange(desc(AUC))
# res = data.frame(res)
# 
# 
# benchmark$Approach = factor(benchmark$Approach,levels = res$Approach)
# res = benchmark %>% 
#   group_by(Approach,Dataset,Seed) %>% 
#   summarise_all(mean) %>% 
#   arrange(desc(AUC))
# res = data.frame(res)
# write_csv(benchmark,paste0(f_name,".csv"))
# 
# tiff(filename =paste0(f_name,".tiff"),width = 4.5,height = 5.5,units = "in",res = 300)
# ggplot(benchmark,aes(Approach,AUC))+
#   theme_bw()+
#   coord_flip()+
#   stat_summary(fun.y = mean, geom = "point",size = 5,col = "black")+
#   stat_summary(fun.data = mean_se,geom = "errorbar")+
#   theme(legend.position = "top",
#         plot.title = element_text(size = 7,hjust = .5,face = "bold"),
#         #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
#         axis.title = element_text(size = 12),
#         #axis.title.y = element_blank(),
#         #axis.text.y = element_text(size = 7),
#         #axis.text.y = element_blank(),
#         #legend.margin=margin(-1,-1,-1,-1),
#         strip.switch.pad.wrap = margin(0,0,0,0),
#         legend.margin=margin(-5,-10,-10,-10),
#         axis.text = element_text(size = 12),
#         #panel.grid = element_blank(),
#         legend.key.size = unit(.15,units = "in"),
#         legend.text = element_text(size = 8),
#         legend.title = element_text(size = 8),
#         #legend.background = element_rect(colour = "black")
#   )
# 
# dev.off()
# 
# 
# 
# 
# ## kruskal test
# kw = spread(benchmark[,1:6],"Approach","AUC")
# kruskal.test(x = res$AUC,g = res$Approach)
# 
# wilcox.test(x = kw$`DCV-rfRFE`,y = kw$`CLR-LASSO`,paired = T,alternative = "two.sided")
# wilcox.test(x = kw$`DCV-ridgeEnsemble`,y = kw$`Coda-LASSO`,paired = T,alternative = "two.sided")
