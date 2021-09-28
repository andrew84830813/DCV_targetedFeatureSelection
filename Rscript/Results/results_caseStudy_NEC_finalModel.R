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






## read data
load("Output/microbiomeDB_dataAtLeast7DaysBeforeNEC.Rda")
df = exp$mbiome
md1 = exp$metadata



## View Class Labels
table(df$Status)
## Get folds
fld = rep(1,nrow(df))

benchmark = data.frame()



td = selEnergyPermR::processCompData(df, minPrevalence = .9)
impFact = td$impFactor
minorityClass = td$minClss
majorityClass = td$majClass
dat = td$processedData
y_labels = dat[,1]
dat = dat[,-1]

## Merge with samples that have mbiome data
md1 = left_join(data.frame(X1 = rownames(dat)),md1)
  
  
  # DCV --------------------------------------------------------------
  
  perc_totalParts2Keep = .75
  num_sets = 5
  
  base_dims = ncol(dat)
  max_parts = round(perc_totalParts2Keep*base_dims)
  sets = round(seq(1,max_parts,length.out = num_sets))[-1]
  
  
  # Run K-Fold Cross Validation ---------------------------------------------
  inner_perf = data.frame()
  
  
  ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")
  #ensemble = c("ranger","xgbTree","xgbLinear")
  max_sparsity = .9
  classes = as.character(unique(y_labels))
  
  
  
  ## Tune target features
  for(sd1 in 1:1){
    set.seed(sd1)
    k_fold = 2
    overll_folds = caret::createFolds(y_labels,k = k_fold,list = F)
    innerfold_data = lodo_partition(data.frame(Status = y_labels,dat),
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
  trainx = data.frame(fastImputeZeroes(dat,impFactor = impFact))
  testx = data.frame(fastImputeZeroes(dat,impFactor = impFact)) 
  
  ## compute log ratios
  lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = y_labels,trainx))
  lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = y_labels,testx))
  
  cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train, 
                                      includeInfoGain = T, nfolds = 1, numRepeats = 1, 
                                      rankOrder = F)
  
  ## Apply targted feature selection method
  tar_Features = inner_perf2$tar_Features[which.max(inner_perf2$AUC)]
  
  tar_dcv = targeted_dcvSelection(trainx = trainx,useRidgeWeights = T,scaledata = F,
                                  testx = testx,
                                  dcv = cc.dcv$lrs,lrs.train = lrs.train,lrs.test = lrs.test,
                                  y_label = y_labels,
                                  seed = sd,
                                  ensemble = ensemble,
                                  y_test = y_labels,
                                  tarFeatures = tar_Features,
                                  ts.id = data.frame(ID= 1:nrow(dat),Status = y_labels), 
                                  max_sparsity = max_sparsity
  )
  
  
## Save Results
#save(tar_dcv,file = "Output/caseStudy_NEC_finalModel.Rda")


##Load Results
load("Output/caseStudy_NEC_finalModel.Rda")
library(ggsignif)

##Retrieve Coeff
cof = coef(tar_dcv$glm_model$mdl)
cc = tar_dcv$ridge_coefficients
p = stats::predict(tar_dcv$glm_model$mdl, newx = as.matrix(tar_dcv$glm_model$data$train), s = "lambda.min",type = "lin")
p.df = data.frame(NEC = y_labels,GLM_Score = as.numeric(p)-cof[1],
                  Survival = factor(md1$Survived,levels = c("Yes","No")),
                  Sepsis = md1$Sepsis.diagnosed,onsetDays = md1$Days.of.period.NEC.diagnosed)


## Visualize Overall Scores Between Groups
pdf(file = "Figures/caseStudy_NEC_glmScores.pdf",width = 2 ,height = 2.5)
ggplot(p.df,aes(NEC,GLM_Score,fill = NEC))+
  geom_boxplot(alpha = .7)+
  #geom_jitter(aes(col  = NEC),width = .1,alpha = .9)+
  ggsci::scale_fill_lancet()+
  ggsci::scale_color_lancet()+
  ylab("Logistics Regression Score")+
  # geom_signif(comparisons = list(c("No", "Yes")), 
  #             map_signif_level=F,test = "wilcox.test",)+
  theme_bw()+
  theme(legend.position = "none",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 8,face = "bold"),
        #axis.title.y = element_blank(),
        #axis.title.x = element_text(size = 8,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+
        
  )
dev.off()

## Relate Scores to Clinical Outcomes
p.df1 = p.df %>% 
  filter(NEC=="Yes")

pdf(file = "Figures/caseStudy_NEC_clinOutcome_Surv.pdf",width = 2.2 ,height = 2.25)
ggplot(p.df1,aes(Survival,GLM_Score))+
  geom_boxplot(alpha = .7)+
  #geom_jitter(aes(col  = NEC),width = .1,alpha = .9)+
  ggsci::scale_fill_lancet()+
  ggsci::scale_color_lancet()+
  ylab("Logistics Regression Score")+
  geom_signif(comparisons = list(c("No", "Yes")),tip_length = 0, 
              map_signif_level=F,test = "wilcox.test",)+
  theme_bw()+
  theme(legend.position = "none",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 8,face = "bold"),
        #axis.title.y = element_blank(),
        #axis.title.x = element_text(size = 8,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+
        
  )
dev.off()

## Relate Scores to time
p.df2 = data.frame(p.df1,bin = cut(p.df1$onsetDays,breaks = seq(-7,-25,by = -1) ))
table(p.df2$bin)
levels(p.df2$bin) = rev(levels(p.df2$bin))
p.df2$bin = factor(p.df2$bin,labels = 7:24)
p.df2 = p.df2 %>% 
  filter(!is.na(bin))
p.df1.null = p.df %>% 
  filter(NEC!="Yes")

pdf(file = "Figures/caseStudy_NEC_clinOutcome_onset.pdf",width = 3.25 ,height = 2.5)
ggplot(p.df2,aes(bin,GLM_Score))+
  #geom_point(alpha = .7)+
  stat_summary(fun.y = mean, geom = "line",col = "black",aes(group =1))+
  stat_summary(fun.y = mean, geom = "point",size = 2,col = "black")+
  stat_summary(fun.data = mean_se,geom = "errorbar",width = .35)+
  geom_hline(yintercept = mean(p.df1.null$GLM_Score),lty = "dashed",col = "red")+
  #geom_jitter(aes(col  = NEC),width = .1,alpha = .9)+
  ggsci::scale_fill_lancet()+
  ggsci::scale_color_lancet()+
  xlab("NEC Onset (Days)")+
  ylab("Logistics Regression Score")+
  geom_signif(comparisons = list(c("No", "Yes")),tip_length = 0, 
              map_signif_level=F,test = "wilcox.test",)+
  theme_bw()+
  theme(legend.position = "none",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 8,face = "bold"),
        #axis.title.y = element_blank(),
        #axis.title.x = element_text(size = 8,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+
        
  )

dev.off()


## Process Weights and construct log ratio network
library(igraph)
imp.df = data.frame(Ratio = names(cc),Imp = abs(as.numeric(cc)),raw = as.numeric(cc))
keyRats = tidyr::separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)
## Define weight such that:  weight * log(a/b) = weight * log(a) - weight * log(b) 
weights.df = data.frame(Part = keyRats$Num,Coef = keyRats$raw)
weights.df = rbind(weights.df,data.frame(Part = keyRats$Denom,Coef = -1*keyRats$raw))
weights.df = weights.df %>% 
  group_by(Part) %>% 
  summarise_all(.funs = sum)
weights.df$col = if_else(weights.df$Coef>0,"NEC","nonNEC")
weights.df$col = factor(weights.df$col,levels = c("nonNEC","NEC"))

library(igraph)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
el_act = data.frame(get.edgelist(g))
el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
el_act = left_join(el_act,imp.df)


vertices = data.frame(Part = V(g)$name,Label =  V(g)$name)
vertices = left_join(vertices,weights.df)
v_color = if_else(vertices$Coef>0,ggsci::pal_lancet(alpha = .9)(2)[1],ggsci::pal_lancet(alpha = .9)(2)[2])
vertices$abs_coef = abs(vertices$Coef)

el_act = data.frame(get.edgelist(g))
el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
el_act = left_join(el_act,imp.df)
col = if_else(el_act$raw>0,ggsci::pal_lancet(alpha = .2)(2)[1],ggsci::pal_lancet(alpha = .2)(2)[2])
el_act$col = col


E(g)$weight = el_act$Imp
toGephi(Graph = g,Name = "NEC_net",attribute_metadata = vertices,edge_metadata = el_act)


## Visualize Net Contribution of Parts to Score
pdf(file = "Figures/caseStudy_NEC_feat_contrib.pdf",width = 5 ,height = 4.5)
ggplot(weights.df,aes(reorder(Part,Coef),Coef,fill = col))+
  geom_col(width = .7)+
  ggsci::scale_fill_aaas()+
  facet_wrap(.~col,scales = "free")+
  coord_flip()+
  #scale_fill_manual(values = ggsci::pal_aaas()(2)[2:1])+
  theme_bw()+
  theme()+
  theme(legend.position = "none",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        strip.background = element_blank(),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        #axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+
        
  )
dev.off()





## Visulaize Log Ratio Network of GLM Coeff 
## Can think of as a microbial network
pdf(file = "Figures/caseStudy_CRC_lrnet.pdf",width = 5 ,height = 5)
par(mar=c(0,0,0,0)+.1)
plot(g,
     vertex.color = v_color,#rep(ggsci::pal_lancet(alpha = .75)(9)[8],vcount(g)),
     layout = igraph::layout.kamada.kawai,
     vertex.frame.color = "white",vertex.frame.width = .5,
     vertex.label.cex = .7,
     vertex.label.color = "black",
     edge.color  =col,
     edge.width = 125*el_act$Imp ,
     vertex.size = abs(vertices$Coef) * 15+2,
     edge.curved = .2,
     edge.arrow.size = .51)
dev.off()








































## Process Weights and construct log ratio network
library(igraph)
imp.df = data.frame(Ratio = names(cc),Imp = abs(as.numeric(cc)),raw = as.numeric(cc))
keyRats = tidyr::separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)
## Define weight such that:  weight * log(a/b) = weight * log(a) - weight * log(b) 
weights.df = data.frame(Part = keyRats$Num,Coef = keyRats$raw)
weights.df = rbind(weights.df,data.frame(Part = keyRats$Denom,Coef = -1*keyRats$raw))
weights.df = weights.df %>% 
  group_by(Part) %>% 
  summarise_all(.funs = sum)
weights.df$col = if_else(weights.df$Coef>0,"NEC","Control")
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
# g = igraph::simplify(g, remove.loops = TRUE,
#                      edge.attr.comb = igraph_opt("edge.attr.comb"))
el_act = data.frame(get.edgelist(g))
el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
el_act = left_join(el_act,imp.df)
col = if_else(el_act$raw>0,ggsci::pal_lancet(alpha = .7)(2)[1],ggsci::pal_lancet(alpha = .7)(2)[2])

vertices = data.frame(Part = V(g)$name)
vertices = left_join(vertices,weights.df)
v_color = if_else(vertices$Coef>0,ggsci::pal_lancet(alpha = .5)(2)[2],ggsci::pal_lancet(alpha = .5)(2)[1])



## Visualize Net Contribution of Parts to Score
ggplot(weights.df,aes(reorder(Part,Coef),Coef,fill = col))+
  geom_col(col = "black")+
  ggsci::scale_fill_aaas()+
  coord_flip()+
  theme_bw()+
  theme()+
  theme(legend.position = "top",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        #axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+
        
  )



## Visulaize Log Ratio Network of GLM Coeff 
## Can think of as a microbial network
pdf(file = "Figures/caseStudy_NEC_lrnet.pdf",width = 5 ,height = 5)
par(mar=c(0,0,0,0)+.1)
plot(g,
     vertex.color = v_color,#rep(ggsci::pal_lancet(alpha = .75)(9)[8],vcount(g)),
     layout = igraph::layout.kamada.kawai,
     vertex.frame.color = "white",
     vertex.label.cex = .5,
     vertex.label.color = "black",
     edge.color  =col,
     edge.width = 125*el_act$Imp ,
     vertex.size = abs(vertices$Coef) *150+5,
     #edge.curved = .2,
     edge.arrow.size = 1)
dev.off()





















































# 
# 
# 
# 
# 
# 
# 
# 
# 
# ### Metada Data
# ## Pre-Process Metadata
# mm.char = md1 %>% 
#   select(Sex,
#          Chorioamnionitis,
#          Maternal.antepartum.antibiotics.administered,
#          Baby.s.delivery.mode,
#          Born.stunted,
#          Host.diet)
# ##-----------------------------
# ## Add covairate
# ## cat features
# dummy <- dummyVars(" ~ .", data=rbind(mm.char))
# newdata <- data.frame(predict(dummy, newdata = rbind(mm.char)))
# newdata.train = newdata[1:nrow(mm.char),]
# 
# 
# ## cont features
# train_metadata = data.frame(Age = md1$Age.at.sample.collection..days.,
#                             gesAge = md1$Gestational.age.at.birth..days.,
#                             host_weight = md1$Host.weight..g.,
#                             newdata.train
# )
# 
# 
# 
# # GLM - Meta --------------------------------------------------------------------
# message("Compute Performance - Metadata Alone GLM")
# suppressMessages(suppressWarnings({
#   
#   
#   ## retrieve test and train data
#   train.data = cbind(train_metadata)
#   test.data = cbind(train_metadata)
#   y_label = y_labels
#   y_test = y_labels
#   
#   
#   ## Apply Penalized Regression
#   ## Tune Alpha
#   type_family = if_else(length(classes)>2,"multinomial","binomial")
#   infld = 2
#   flds = caret::createFolds(y = y_label,k = infld,list = F)
#   compTime2 = system.time({
#     aseq = seq(1e-3,1,length.out = 10)
#     min_dev =  foreach(a = aseq,.combine = rbind)%dopar%{
#       
#       aucc = c()
#       for(f in 1:infld){
#         bool  = flds==f
#         compTime2 = system.time({
#           cv.clrlasso <- glmnet::cv.glmnet(as.matrix(train.data[bool,]),y_label[bool], 
#                                            standardize=T, alpha=a,family=type_family)
#         })
#         
#         ## make predictions
#         p = predict(cv.clrlasso, newx = as.matrix(train.data[!bool,]), s = "lambda.min",type = "response")
#         if(type_family=="binomial"){
#           mroc = pROC::roc(y_label[!bool],p)
#           mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
#         }else{
#           ## multiclass
#           mroc = pROC::multiclass.roc(y_label[!bool],p[,,1])
#           mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
#         }
#         aucc = c(aucc,as.numeric(mroc.dcvlasso))
#       }
#       data.frame(a,auc = mean(aucc))
#     }
#   })
#   min_dev = min_dev %>% 
#     arrange(desc(auc))
#   
#   ## Train GLM
#   compTime2 = system.time({
#     cv.clrlasso <- glmnet::cv.glmnet(as.matrix(train.data),y_label, standardize=T, alpha = min_dev$a[1],family=type_family)
#   })
#   if(type_family=="binomial"){
#     features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
#     features = features[-1,]
#     features = features[abs(features)>0]
#     length(features)
#     c = as.matrix(coef(cv.clrlasso, s = "lambda.min"))[-1,]
#     train_data.metaGLM = subset(train.data,select = names(features))
#     test_data.metaGLM = subset(test.data,select = names(features))
#   }else{
#     features = as.matrix(stats::coef(cv.clrlasso, s = "lambda.min"))
#     feat.df = data.frame()
#     for(o in 1:length(features)){
#       ph = as.matrix(features[[o]])
#       feat = ph[-1,]
#       keep = feat[abs(feat)>0]
#       feat.df = rbind(feat.df,data.frame(Ratio = names(keep),coef = as.numeric(keep)))
#     }
#     feat.df =feat.df %>% 
#       group_by(Ratio) %>% 
#       summarise(coef = sum(coef)) %>% 
#       filter(coef!=0)
#     train_data.metaGLM = subset(train.data,select = feat.df$Ratio)
#     test_data.metaGLM = subset(test.data,select = feat.df$Ratio)
#   }
#   
#   
# }))
# 
# 
# 
# # GLM - Meta --------------------------------------------------------------------
# message("Compute Performance - Metadata Alone GLM")
# suppressMessages(suppressWarnings({
#   
#   
#   ## retrieve test and train data
#   cc = tar_dcv$ridge_coefficients
#   train.data = cbind(sweep(tar_dcv$weighted_features$train,MARGIN = 2,STATS = as.numeric(cc),FUN = "/")  ,train_data.metaGLM)
#   test.data = cbind(sweep(tar_dcv$weighted_features$test,MARGIN = 2,STATS = as.numeric(cc),FUN = "/"),test_data.metaGLM)
#   y_label = y_labels
#   y_test = y_labels
#   
#   
#   ## Apply Penalized Regression
#   ## Tune Alpha
#   type_family = if_else(length(classes)>2,"multinomial","binomial")
#   infld = 2
#   flds = caret::createFolds(y = y_label,k = infld,list = F)
#   compTime2 = system.time({
#     aseq = seq(1e-3,1,length.out = 10)
#     min_dev =  foreach(a = aseq,.combine = rbind)%dopar%{
#       
#       aucc = c()
#       for(f in 1:infld){
#         bool  = flds==f
#         compTime2 = system.time({
#           cv.clrlasso <- glmnet::cv.glmnet(as.matrix(train.data[bool,]),y_label[bool], 
#                                            standardize=T, alpha=a,family=type_family)
#         })
#         
#         ## make predictions
#         p = predict(cv.clrlasso, newx = as.matrix(train.data[!bool,]), s = "lambda.min",type = "response")
#         if(type_family=="binomial"){
#           mroc = pROC::roc(y_label[!bool],p)
#           mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
#         }else{
#           ## multiclass
#           mroc = pROC::multiclass.roc(y_label[!bool],p[,,1])
#           mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
#         }
#         aucc = c(aucc,as.numeric(mroc.dcvlasso))
#       }
#       data.frame(a,auc = mean(aucc))
#     }
#   })
#   min_dev = min_dev %>% 
#     arrange(desc(auc))
#   
#   ## Train GLM
#   compTime2 = system.time({
#     cv.clrlasso <- glmnet::cv.glmnet(as.matrix(train.data),y_label, standardize=T, alpha = min_dev$a[1],family=type_family)
#   })
#   if(type_family=="binomial"){
#     features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
#     features = features[-1,]
#     features = features[abs(features)>0]
#     length(features)
#     c = as.matrix(coef(cv.clrlasso, s = "lambda.min"))[-1,]
#   }else{
#     features = as.matrix(stats::coef(cv.clrlasso, s = "lambda.min"))
#     feat.df = data.frame()
#     for(o in 1:length(features)){
#       ph = as.matrix(features[[o]])
#       feat = ph[-1,]
#       keep = feat[abs(feat)>0]
#       feat.df = rbind(feat.df,data.frame(Ratio = names(keep),coef = as.numeric(keep)))
#     }
#     feat.df =feat.df %>% 
#       group_by(Ratio) %>% 
#       summarise(coef = sum(coef)) %>% 
#       filter(coef!=0)
#   }
#   
#   ## make predictions
#   p = predict(cv.clrlasso, newx = as.matrix(test.data), s = "lambda.min",type = "response")
#   if(type_family=="binomial"){
#     mroc = pROC::roc(y_test,p)
#     mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
#   }else{
#     ## multiclass
#     mroc = pROC::multiclass.roc(y_test,p[,,1])
#     mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
#   }
#   
#   ## Compute Number of Part
#   cn =names(c[abs(c)>0])
#   n_ratios = length(cn)
#   uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
#   n_parts  = dplyr::n_distinct(uniqueParts)
#   
#   ## Save Performance
#   perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
#                     Dataset = f_name,Seed = sd,Fold = f,Approach = "meta_dataGLM",
#                     AUC = as.numeric(pROC::auc(mroc)),
#                     number_parts = n_parts,number_ratios = ncol(train.data) ,comp_time = NA,
#                     base_dims = ncol(train.data))
#   benchmark = rbind(benchmark,perf)
#   
#   
# }))
# 
# 
# 
# ## get coeff
# ## Accounting for metadata
# c = as.matrix(coef(cv.clrlasso, s = "lambda.min"))[-1,]
# 
# cf = c[abs(c)>0]
# cf = cf[str_detect(string = names(cf),pattern = "___")]
# imp.df = data.frame(Ratio = names(cf),Imp = abs(as.numeric(cf)))
# keyRats = tidyr::separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)
# weights.df = data.frame(Part = keyRats$Num,Coef = keyRats$Imp)
# weights.df = rbind(weights.df,data.frame(Part = keyRats$Denom,Coef = cf))
# weights.df = weights.df %>% 
#   group_by(Part) %>% 
#   summarise_all(.funs = sum)
# weights.df$col = if_else(weights.df$Coef>0,"Control","NEC")
# ggplot(weights.df,aes(reorder(Part,Coef),Coef,fill = col))+
#   geom_col()+
#   ggsci::scale_fill_aaas()+
#   coord_flip()+
#   theme_bw()
# 
# 
# 
# 
# p = predict(cv.clrlasso, newx = as.matrix(test.data), s = "lambda.min",type = "link")
# p.df = data.frame(Status = y_labels,Res = as.numeric(p))
# ggplot(p.df,aes(Status,Res,fill  = Status))+
#   geom_boxplot(alpha = .6)+
#   geom_jitter(aes(col  =Status),width = .1,alpha = .9)+
#   ggsci::scale_fill_aaas()+
#   ggsci::scale_color_aaas()
# 
# ggplot(weights.df,aes(reorder(Part,Coef),Coef,fill = col))+
#   geom_col()+
#   ggsci::scale_fill_aaas()+
#   coord_flip()+
#   theme_bw()
# 
# 
