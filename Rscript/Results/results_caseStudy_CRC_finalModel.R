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






## load data
load("Output/LODO/crc_data.Rda")
df = obj$data
md = obj$md





## View Class Labels
table(df$Status)
combs = data.frame(combs)

## get read of samples with reads <10e6
rs = rowSums(df[,-1])
bool = rs>10e6
df = df[bool,]
md = md[bool,]
table(md$dataset)
benchmark = data.frame()




td = selEnergyPermR::processCompData(df, minPrevalence = .9)
impFact = td$impFactor
minorityClass = td$minClss
majorityClass = td$majClass
dat = td$processedData
y_labels = dat[,1]
dat = dat[,-1]



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
tar_Features = 50

tar_dcv = targeted_dcvSelection(trainx = trainx,
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
#save(tar_dcv,file = "Output/caseStudy_CRC_finalModel_top50.Rda")


##Load Results
load("Output/caseStudy_CRC_finalModelopt.Rda")
load("Output/caseStudy_CRC_finalModel_top50.Rda")

library(ggsignif)

##Retrieve Coeff
cof = coef(tar_dcv$glm_model$mdl)
cc = tar_dcv$ridge_coefficients
sample_id = rownames(df)
p = stats::predict(tar_dcv$glm_model$mdl, newx = as.matrix(tar_dcv$glm_model$data$train), s = "lambda.min",type = "lin")
p.df = data.frame(CRC = y_labels,GLM_Score = as.numeric(p),Sample_ID = sample_id)


## Visualize Overall Scores Between Groups
pdf(file = "Figures/caseStudy_CRC_scores.pdf",width = 2.5 ,height = 2.2)
ggplot(p.df,aes(CRC,GLM_Score-cof[1],fill = CRC))+
  geom_boxplot(alpha = .7,outlier.color = NA)+
  #geom_jitter(aes(col  = CRC),width = .05,alpha = .5)+
  ggsci::scale_fill_lancet()+
  ggsci::scale_color_lancet()+
  ylab("Logistics Regression Score")+
  # geom_signif(comparisons = list(c("control", "CRC")),
  #             map_signif_level=F,test = "wilcox.test")+
  theme_bw()+
  theme(legend.position = "none",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 8,face = "bold"),
        #axis.title.y = element_blank(),
        axis.title.x = element_blank(),
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
# AJCC Staging Yachida Datset ---------------------------------------------

load("Output/CRC_exp/curatedMetaGenome_YachidaS_2019.Rda")
clin_md = expResults$metadata
clin_md = inner_join(data.frame(data.frame(sampleID = sample_id,CRC = y_labels,GLM_Score = as.numeric(p))),clin_md)
clin_md = clin_md %>% 
  filter(!is.na(ajcc))

# ggplot(clin_md,aes(disease_subtype,GLM_Score-cof[1]))+
#   geom_boxplot(alpha = .7)+
#   #geom_jitter(aes(col  = NEC),width = .1,alpha = .9)+
#   ggsci::scale_fill_lancet()+
#   ggsci::scale_color_lancet()+
#   ylab("Score")+
#   # geom_signif(comparisons = list(c("No", "Yes")),tip_length = 0, 
#   #             map_signif_level=F,test = "wilcox.test",)+
#   theme_bw()+
#   theme(legend.position = "none",panel.grid = element_blank(),
#         plot.title = element_text(size = 8,hjust = .5,face = "bold"),
#         #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
#         axis.title = element_text(size = 8,face = "bold"),
#         #axis.title.y = element_blank(),
#         #axis.title.x = element_text(size = 8,face = "bold"),
#         #axis.text.y = element_text(size = 7),
#         #axis.text.y = element_blank(),
#         #legend.margin=margin(-1,-1,-1,-1),
#         strip.switch.pad.wrap = margin(0,0,0,0),
#         legend.margin=margin(-5,-10,-10,-10),
#         axis.text = element_text(size = 8),
#         #panel.grid = element_blank(),
#         legend.key.size = unit(.15,units = "in"),
#         legend.text = element_text(size = 8),
#         legend.title = element_blank(),
#         #legend.background = element_rect(colour = "black")+
#         
#   )

stage = data.frame(ajcc = unique(clin_md$ajcc),cStage = c("SIII/IV","SI/II","SIII/IV","SI/II","S0"))
clin_md = left_join(clin_md,stage)
clin_md$cStage = factor(clin_md$cStage)
clin_md$Score = clin_md$GLM_Score - cof[1]
table(clin_md$cStage)

## Linear Model
md = lm(Score~cStage,data = clin_md)
sm  = (anova(md))

pdf(file = "Figures/caseStudy_CRC_clinOutcome.pdf",width = 2.5 ,height = 2.2)
ggplot(clin_md,aes(cStage,GLM_Score-cof[1]))+
  geom_boxplot(alpha = .7)+
  #geom_jitter(aes(col  = NEC),width = .1,alpha = .9)+
  ggsci::scale_fill_lancet()+
  ggsci::scale_color_lancet()+
  xlab("AJCC Stage")+
  ylab("Logistics Regression Score")+
  labs(caption = paste0("F=",round(sm$`F value`[1],5),", p = ",round(sm$`Pr(>F)`[1],16)))+
  # geom_signif(comparisons = list(c("No", "Yes")),tip_length = 0, 
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





# TNM Staging 6 dtasets -------------------------------------------------------------

stage_data = data.frame()
load("Output/CRC_exp/curatedMetaGenome_YuJ_2015.Rda")
clin_md = expResults$metadata
mt = data.frame(matrix(nrow = 1,ncol = nrow(clin_md),dimnames = list(row = 1,col = clin_md$sampleID)))
clin_md$sampleID = colnames(mt)
clin_md = inner_join(data.frame(data.frame(sampleID = sample_id,CRC = y_labels,GLM_Score = as.numeric(p)-cof[1])),clin_md)
clin_md = clin_md %>% 
  filter(!is.na(tnm)) 
stage = stringr::str_extract((clin_md$tnm), "^.{2}")
ph = data.frame(sampleID = clin_md$sampleID,stage,Score = clin_md$GLM_Score)
stage_data = rbind(stage_data,ph)


load("Output/CRC_exp/curatedMetaGenome_FengQ_2015.Rda")
clin_md = expResults$metadata
mt = data.frame(matrix(nrow = 1,ncol = nrow(clin_md),dimnames = list(row = 1,col = clin_md$sampleID)))
colnames(clin_md)
clin_md$sampleID = colnames(mt)
clin_md = inner_join(data.frame(data.frame(sampleID = sample_id,CRC = y_labels,GLM_Score = as.numeric(p)-cof[1])),clin_md)
clin_md = clin_md %>% 
  filter(!is.na(tnm)) 
stage = stringr::str_extract((clin_md$tnm), "^.{2}")
ph = data.frame(sampleID = clin_md$sampleID,stage,Score = clin_md$GLM_Score)
stage_data = rbind(stage_data,ph)

load("Output/CRC_exp/curatedMetaGenome_ZellerG_2014.Rda")
clin_md = expResults$metadata
mt = data.frame(matrix(nrow = 1,ncol = nrow(clin_md),dimnames = list(row = 1,col = clin_md$sampleID)))
colnames(clin_md)
clin_md$sampleID = colnames(mt)
clin_md = inner_join(data.frame(data.frame(sampleID = sample_id,CRC = y_labels,GLM_Score = as.numeric(p)-cof[1])),clin_md)
clin_md = clin_md %>% 
  filter(!is.na(tnm)) 
stage = stringr::str_extract((clin_md$tnm), "^.{2}")
ph = data.frame(sampleID = clin_md$sampleID,stage,Score = clin_md$GLM_Score)
stage_data = rbind(stage_data,ph)

# load("Output/CRC_exp/curatedMetaGenome_ThomasAM_2019_c.Rda")
clin_md = expResults$metadata
mt = data.frame(matrix(nrow = 1,ncol = nrow(clin_md),dimnames = list(row = 1,col = clin_md$sampleID)))
colnames(clin_md)
clin_md$sampleID = colnames(mt)
clin_md = inner_join(data.frame(data.frame(sampleID = sample_id,CRC = y_labels,GLM_Score = as.numeric(p)-cof[1])),clin_md)
clin_md = clin_md %>%
  filter(!is.na(tnm))
stage = stringr::str_extract((clin_md$tnm), "^.{2}")

ph = data.frame(sampleID = clin_md$sampleID,stage ,Score = clin_md$GLM_Score)
stage_data = rbind(stage_data,ph)

load("Output/CRC_exp/curatedMetaGenome_GuptaA_2019.Rda")
clin_md = expResults$metadata
mt = data.frame(matrix(nrow = 1,ncol = nrow(clin_md),dimnames = list(row = 1,col = clin_md$sampleID)))
colnames(clin_md)
clin_md$sampleID = colnames(mt)
clin_md = inner_join(data.frame(data.frame(sampleID = sample_id,CRC = y_labels,GLM_Score = as.numeric(p)-cof[1])),clin_md)
clin_md = clin_md %>% 
  filter(!is.na(tnm)) 
stage = stringr::str_extract((clin_md$tnm), "^.{2}")
ph = data.frame(sampleID = clin_md$sampleID,stage,Score = clin_md$GLM_Score)
stage_data = rbind(stage_data,ph)


load("Output/CRC_exp/curatedMetaGenome_WirbelJ_2018.Rda")
clin_md = expResults$metadata
mt = data.frame(matrix(nrow = 1,ncol = nrow(clin_md),dimnames = list(row = 1,col = clin_md$sampleID)))
colnames(clin_md)
clin_md$sampleID = colnames(mt)
clin_md = inner_join(data.frame(data.frame(sampleID = sample_id,CRC = y_labels,GLM_Score = as.numeric(p)-cof[1])),clin_md)
clin_md = clin_md %>% 
  filter(!is.na(tnm)) 
stage = stringr::str_extract((clin_md$tnm), "^.{2}")
ph = data.frame(sampleID = clin_md$sampleID,stage,Score = clin_md$GLM_Score)
stage_data = rbind(stage_data,ph)

#stage = data.frame(tnm = unique(clin_md$tnm),Stage = stage)
stage = data.frame(tnm = unique(clin_md$tnm),Stage = if_else(stage%in%c("t1","t2"),"t1/2","t3/4"))
clin_md = left_join(clin_md,stage)
clin_md$Score = clin_md$GLM_Score - cof[1]


stage_data = stage_data %>% 
  filter(stage %in% c("t1","t2","t3","t4"))


## Linear Model
table(stage_data$stage)
md = lm(Score~stage,data = stage_data)
sm = (anova(md))

pdf(file = "Figures/caseStudy_CRC_clinOutcomeAllData.pdf",width = 2.5 ,height = 2.2)
ggplot(stage_data,aes(stage,Score))+
  geom_boxplot(alpha = .7)+
  #geom_jitter(aes(col  = NEC),width = .1,alpha = .9)+
  ggsci::scale_fill_lancet()+
  ggsci::scale_color_lancet()+
  xlab("TNM Stage")+
  ylab("Logistics Regression Score")+
  labs(caption = paste0("F=",round(sm$`F value`[1],5),", p = ",round(sm$`Pr(>F)`[1],16)))+
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
weights.df$col = if_else(weights.df$Coef>0,"CRC","Control")


library(igraph)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
el_act = data.frame(get.edgelist(g))
el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
el_act = left_join(el_act,imp.df)



# E(g)$weight = el_act$Imp
# adjacency_matrix <- igraph::as_adjacency_matrix(graph = g, 
#                                                 sparse = F, attr = "weight")
# g = knn_graph(adj_mat = adjacency_matrix,K = 3,plot_TrueFalse = T,sim_ = F)
# g = g$Graph
# E(g)$weight

vertices = data.frame(Part = V(g)$name,Label =  V(g)$name)
vertices = left_join(vertices,weights.df)
v_color = if_else(vertices$Coef>0,ggsci::pal_lancet(alpha = .9)(2)[2],ggsci::pal_lancet(alpha = .9)(2)[1])
vertices$abs_coef = abs(vertices$Coef)

el_act = data.frame(get.edgelist(g))
el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
el_act = left_join(el_act,imp.df)
col = if_else(el_act$raw>0,ggsci::pal_lancet(alpha = .2)(2)[2],ggsci::pal_lancet(alpha = .2)(2)[1])
el_act$col = col


#E(g)$weight = el_act$Imp
#toGephi(Graph = g,Name = "CRC_opt",attribute_metadata = vertices,edge_metadata = el_act)


## Visualize Net Contribution of Parts to Score
pdf(file = "Figures/caseStudy_CRC_feat_contrib.pdf",width = 6 ,height = 7.5)
ggplot(weights.df,aes(reorder(Part,Coef),Coef,fill = col))+
  geom_col(width = .7)+
  ggsci::scale_fill_aaas()+
  facet_wrap(.~col,scales = "free")+
  coord_flip()+
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
     vertex.label.cex = .4,
     vertex.label.color = "black",
     edge.color  =col,
     edge.width = 125*el_act$Imp ,
     vertex.size = abs(vertices$Coef) * 15+2,
     edge.curved = .2,
     edge.arrow.size = .51)
dev.off()




























