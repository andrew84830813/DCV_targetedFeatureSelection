
###-------------------------------------------------------*
#### NEC ####
###-------------------------------------------------------*
metaData = data.frame(readr::read_tsv(file = "Data/NICUNEC.WGS.sample_details.tsv"))
otuData = data.frame(readr::read_tsv(file = "Data/NICUNEC.WGS.taxon_abundance.tsv"))


## NEC samples
metaData1 = metaData %>% 
  filter(Days.of.period.NEC.diagnosed<=-7 ) %>% 
  #filter(Days.of.period.NEC.diagnosed<=-7 & Days.of.period.NEC.diagnosed>=-14 ) %>% 
  filter(Age.at.sample.collection..days.>0) 
# 
# fit <- fitdistrplus::fitdist(metaData1$Age.at.sample.collection..days., "nbinom")
# plot(fit)
# fitdistrplus::gofstat(fit,discrete = T)


metaData2 = metaData %>% 
  filter(is.na(Days.of.period.NEC.diagnosed)) %>% 
  filter(Age.at.sample.collection..days. < max(metaData1$Age.at.sample.collection..days.)) 
hist(metaData2$Age.at.sample.collection..days.)
hist(metaData1$Age.at.sample.collection..days.)

## age matched sample
set.seed(08272008)
metaData3 = metaData2

# metaData3 = sample_n(metaData2,size = nrow(metaData1)*2,replace = F,
#                      weight = dnbinom(metaData2$Age.at.sample.collection..days.,size = fit$estimate[1],mu = fit$estimate[2] ))
# hist(metaData3$Age.at.sample.collection..days.)
# hist(metaData1$Age.at.sample.collection..days.)



md1 = rbind(metaData1,metaData3)
md = data.frame(X1 = md1$X1,Status = md1$Necrotizing.enterocolitis)
table(md$Status)
table(md1$Days.of.period.NEC.diagnosed)
hist(md1$Age.at.sample.collection..days.)

samples =md1 %>% 
  group_by(Subject.ID,Necrotizing.enterocolitis) %>% 
  summarise(n = n())

#define caseIDs
otuData = subset(otuData,select = c("X1",metaData$X1))
## Sep
otuData1 = tidyr::separate(otuData,col = 1,
                           into = c("Kingdom","Phylum","Class","Order",
                                    "Family","Genus","Species","fhfhf"),
                           remove = F,sep = ";")

## retain taxa idenfitfed through family level
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
f_name= "tesr"

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

    
    
   
    # DCV Approaches ----------------------------------------------------------
    
    type = "Dense"
    y_label =  ttData$y_train
    dcv_threshold = .25
    tc = .9
    dcv_thres = 0.05
    bootstrap_reps = 1000
    
    ## Should labels be permuted
    if(permute_labels){
      y_label = sample(y_label)
    }
    
    ## dcv lasso
    compTime = system.time({
      ## get training data
      xt =train.data
      xt = data.frame(fastImputeZeroes(xt,impFactor = ttData$imputeFactor))
      xtest = data.frame(fastImputeZeroes(test.data,impFactor = ttData$imputeFactor)) 
      
      
      
      ## compute log ratios
      lrs =calcLogRatio(data.frame(Status =y_label,xt))
      lrs.test =calcLogRatio(data.frame(Status = ttData$y_test,xtest))
      
      
      ## compute dcvSCores
      cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs, 
                                          includeInfoGain = T, 
                                          nfolds = 1, 
                                          numRepeats = 1, 
                                          rankOrder = T
                                          )
      
      ## hypoth
      dcv =cc.dcv$lrs
      w = min(which(dcv$nDistinct==ncol(xt)))
      dcv = dcv[1:w,]
      dcv_str = dcv_strength(dcvScores = dcv)
      ## bootstrap 5 percentile
      bs = foreach::foreach(i =1:bootstrap_reps,.combine = c)%dopar%{
        set.seed(i)
        quantile(sample(dcv_str$meanDegree,replace = T),probs = dcv_thres)
      }
      rem1 = dcv_str$Node[dcv_str$meanDegree<=median(bs)]
      rem = matrix(rep(0,length(rem1)*w),nrow = length(rem1))
      for(ll in 1:nrow(rem)){
        rem[ll,] = str_detect(string = dcv$Ratio,pattern = rem1[ll])
      }
      bool = if_else(colSums(rem)>0,F,T)
      keep = dcv$Ratio[bool]
      glm.train = getLogratioFromList(keep,raMatrix = xt,Class = "train")
      glm.test = getLogratioFromList(keep,raMatrix = xtest,Class = "test")
      
      
      
    })
    
    # 
    # mm = sapply(train.md, typeof)
    # mm.char = which(mm=="character")[-1]
    # mm.char = train.md[,mm.char]
    # keep= sapply(1:ncol(mm.char), function(x) sum(is.na(mm.char[,x])))
    # mm.char = mm.char[,which(keep==0)]
    # keep= sapply(1:ncol(mm.char), function(x) n_distinct(mm.char[,x])  )
    # mm.char = mm.char[,which(keep>1)]
    
    #colnames(mm.char)
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
    
    # for(i in 1:ncol(newdata.train)){
    #   newdata.train[,i] = factor(newdata.train[,i])
    #   newdata.test[,i] = factor(newdata.test[,i])
    # }

    ## cont features
    train_metadata = data.frame(Age = train.md$Age.at.sample.collection..days.,
                               gesAge = train.md$Gestational.age.at.birth..days.,
                               host_weight = train.md$Host.weight..g.,
                               newdata.train
    )
    glm.train1 = glm.train
    glm.train = cbind(glm.train,train_metadata)
    
    
    test.metadata = data.frame(Age = test.md$Age.at.sample.collection..days.,
                               gesAge = test.md$Gestational.age.at.birth..days.,
                               host_weight = test.md$Host.weight..g.,
                               newdata.test
    )
    glm.test1 = glm.test
    glm.test = cbind(glm.test,test.metadata)
    y_label1 = y_label
    y_label2 = y_label
    
    ##---------------------------
    # 
    # # ## apply smote
    # imbal_train = data.frame(Status = y_label,data.frame((glm.train)))
    # imbal_train <- ROSE::ROSE(Status ~ ., data  = imbal_train,)$data
    # glm.train = data.frame((imbal_train[,-1]))
    # y_label = imbal_train[,1]
    # 
    # ## mbiome alone
    # imbal_train = data.frame(Status = ttData$y_train,data.frame((glm.train1)))
    # imbal_train <- ROSE::ROSE(Status ~ ., data  = imbal_train,)$data
    # glm.train1 = data.frame((imbal_train[,-1]))
    # y_label1 = imbal_train[,1]
    # 
    # ## metadata alone
    # imbal_train = data.frame(Status = ttData$y_train,data.frame((train_metadata)))
    # imbal_train <- ROSE::ROSE(Status ~ ., data  = imbal_train,)$data
    # train_metadata = data.frame((imbal_train[,-1]))
    # y_label2 = imbal_train[,1]
    # 

    
    
    ##scale data
    pp = caret::preProcess(glm.train,method = "scale")
    glm.train <- predict(pp, glm.train)
    glm.test     <- predict(pp, glm.test)
    
    
    pp = caret::preProcess(glm.train1,method = "scale")
    glm.train1 <- predict(pp, glm.train1)
    glm.test1     <- predict(pp, glm.test1)
    
    
    
  
    
    
    # RFE  --------------------------------------------------------------------
    ensemble = c("ranger","xgbTree","gbm")
    suppressMessages(suppressWarnings({
      
      # type = "Dense"
      # ## Feature Selection: rf-RFE 
      # tc = .99
      # c1.cor = cor(glm.train,method = "spearman")
      # c.fc = data.frame(Ratio = caret::findCorrelation(c1.cor,cutoff = tc,names = T))
      # keep =!colnames(glm.train)%in%c.fc$Ratio
      # glm.train = glm.train[,keep]
      # glm.test = glm.test[,keep]
      
      compTime2 = system.time({
        pp = rfeSelection.ByMetric(train_ratio = train_metadata,
                                   test_ratio = test.metadata,
                                   ytrain =y_label2,
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
      #pmat[,4:ncol(pmat)] = clo((pmat[,4:ncol(pmat)]))
      mroc = pROC::multiclass.roc(pmat$Status,pmat[,classes]);mroc
      
      
      ## Compute Number of Part
      cn = colnames(train_data.meta)
      n_ratios = length(cn)
      uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
      n_parts  = dplyr::n_distinct(uniqueParts)
      
      ## Save Performance
      perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                        Dataset = f_name,Seed = sd,Fold = f,Approach = "meta_alone",AUC = as.numeric(pROC::auc(mroc)),
                        number_parts = n_parts,number_ratios = ncol(train_data.meta) ,comp_time = compTime[3]+compTime2[3],
                        base_dims = ncol(train.data))
      benchmark = rbind(benchmark,perf)
    }))
    
    
    
    # RFE  --------------------------------------------------------------------
    suppressMessages(suppressWarnings({
      
      # type = "Dense"
      # ## Feature Selection: rf-RFE 
      # tc = .99
      # c1.cor = cor(glm.train,method = "spearman")
      # c.fc = data.frame(Ratio = caret::findCorrelation(c1.cor,cutoff = tc,names = T))
      # keep =!colnames(glm.train)%in%c.fc$Ratio
      # glm.train = glm.train[,keep]
      # glm.test = glm.test[,keep]
      
      compTime2 = system.time({
        pp = rfeSelection.ByMetric(train_ratio = glm.train1,
                                   test_ratio = glm.test1,
                                   ytrain =y_label1,
                                   ntrees = 750,
                                   sets = 10,
                                   impMeasure = "impurity_corrected",
                                   kfold = 5,
                                   minPercentFeatReturn = .3)
      })
      train_data2 = pp$reducedTrainRatios
      test_data2 = pp$reducedTestRatio
      message("number of features = ",ncol(train_data2))
      keep.mbiome = colnames(train_data2)
      
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
                        Dataset = f_name,Seed = sd,Fold = f,Approach = "mbiome_alone",AUC = as.numeric(pROC::auc(mroc)),
                        number_parts = n_parts,number_ratios = ncol(train_data2) ,comp_time = compTime[3]+compTime2[3],
                        base_dims = ncol(train.data))
      benchmark = rbind(benchmark,perf)
    }))
    
    
    
    # RFE  --------------------------------------------------------------------
    suppressMessages(suppressWarnings({
      
      # type = "Dense"
      # ## Feature Selection: rf-RFE 
      # tc = .99
      # c1.cor = cor(glm.train,method = "spearman")
      # c.fc = data.frame(Ratio = caret::findCorrelation(c1.cor,cutoff = tc,names = T))
      # keep =!colnames(glm.train)%in%c.fc$Ratio
      # glm.train = glm.train[,keep]
      # glm.test = glm.test[,keep]
      
     
      
      train_data2 = subset(glm.train,select = c(keep.mbiome,keep.meta))
      test_data2 = subset(glm.test,select = c(keep.mbiome,keep.meta))
      
      
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
                        number_parts = n_parts,number_ratios = ncol(train_data2) ,comp_time = compTime[3]+compTime2[3],
                        base_dims = ncol(train.data))
      benchmark = rbind(benchmark,perf)
    }))
    
    # RFE  --------------------------------------------------------------------
    suppressMessages(suppressWarnings({
      
      # type = "Dense"
      # ## Feature Selection: rf-RFE 
      # tc = .99
      # c1.cor = cor(glm.train,method = "spearman")
      # c.fc = data.frame(Ratio = caret::findCorrelation(c1.cor,cutoff = tc,names = T))
      # keep =!colnames(glm.train)%in%c.fc$Ratio
      # glm.train = glm.train[,keep]
      # glm.test = glm.test[,keep]
      
  
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
                        Dataset = f_name,Seed = sd,Fold = f,Approach = "mbiome + meta_data.merged",AUC = as.numeric(pROC::auc(mroc)),
                        number_parts = n_parts,number_ratios = ncol(train_data2) ,comp_time = compTime[3]+compTime2[3],
                        base_dims = ncol(train.data))
      benchmark = rbind(benchmark,perf)
    }))
    
    
   
    
    View(benchmark)
    
  }
}

#write_csv(benchmark,file = "Output/case_studyNEC_results.csv")


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
ggplot(res,aes(Approach,AUC))+
  stat_summary(fun.y  = mean,geom = "point")+
  stat_summary(fun.data = mean_se,geom = "errorbar")
  


## kruskal test
kw = spread(benchmark[,1:6],"Approach","AUC")
kruskal.test(x = res$AUC,g = res$Approach)

wilcox.test(x = kw$`mbiome + meta_data`,y = kw$meta_alone,paired = T,alternative = "two.sided")
wilcox.test(x = kw$`DCV-ridgeEnsemble`,y = kw$`Coda-LASSO`,paired = T,alternative = "two.sided")


