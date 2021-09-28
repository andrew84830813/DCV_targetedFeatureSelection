
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
library(ggplot2)
library(ggtern)

## USe to get ggtern to work
#remotes::install_version(package = "ggplot2",version = "3.2.1")


# Load Helper Functions  ---------------------------------------------------------------
fnames = dir("Helper_Functions/")
for(f in fnames){
  source(paste0("Helper_Functions/",f))
}


clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)




## sim parms
cl_center = .5
c2_center = 1
noiseC1 = 2
noiseC2 = -2
nsamps = 100
alr_dimensions = 2

## SIm Dataset
s1 = simDataset(n = nsamps,dims = alr_dimensions,
                offset_c1 = cl_center,
                offset_c2 = c2_center,
                noiseCenter1 = noiseC1,convertCounts = F )
s1 = data.frame(Dataset = "Train",s1)
s2 = simDataset(n = nsamps,dims = alr_dimensions,
                offset_c1 = cl_center,
                offset_c2 = c2_center,
                noiseCenter1 = noiseC2 ,convertCounts = F)
s2 = data.frame(Dataset = "Test",s2)
df = rbind(s1,s2)
tbl = data.frame(df)
tbl$Status = str_replace_all(string = df$Status,pattern = "S",replacement = "C" )



## Overall Data
ggtern(tbl,aes(x = V1,y = V2,z = V3,col = Status, shape = Dataset))+
  geom_point(size = 1,alpha = .5)+
  theme_classic()



## Visulualize Decision Boundary
LR = mean(log(tbl$V1/tbl$V3))
xx = data.frame(Dataset = "",Status = "Expected Decison Bound",V1 = seq(.01,1,length.out = 1000) )
xx$V3 = xx$V1/exp(LR)
xx = xx %>% 
  dplyr::mutate(S = V1+V3) %>% 
  filter(S<=1) %>% 
  dplyr::select(-S) %>% 
  dplyr::mutate(V2 = 1 - (V1 + V3)) %>% 
  dplyr::select(Dataset,Status,V1,V2,V3) 
xxx = rbind(tbl,xx)
ggtern(tbl,aes(x = V1,y = V2,z = V3,col = Status, shape = Dataset))+
  geom_point(size = 2,alpha = .5)+
  geom_line(data = xx,aes(x = V1,y = V2,z = V3))+
  theme_classic()



## Train Model Relative Abundance
## Select Dataset 1
df = tbl %>% 
  filter(Dataset == "Train") %>% 
  select(-Dataset)
df$Status = factor(df$Status)

## Daatset 2
df2 = tbl %>% 
  filter(Dataset == "Test") %>% 
  select(-Dataset)
df2$Status = factor(df2$Status)

## Train RF 
m1 = trainML_Models(trainLRs = df[,-1],testLRs = df2[,-1],
                    ytrain = df[,1],y_test = df2[,1],models = "ranger",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

## Train Performance
m1$performance
## Train Performance
classes = as.character(unique(df$Status))
## Train AUC
preds = m1$models[[1]]$pred
train_mroc = pROC::multiclass.roc(preds$obs,preds[,classes])
train.roc = DiCoVarML::rocPlot(train_mroc)
train.roc$trial = "Train"

## Compute AUC
roc_ = pROC::roc(df2$Status,m1$predictionMatrix[,2])
roc_$thresholds
mroc = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,1:2])
roc.df = DiCoVarML::rocPlot(mroc)
roc.df$trial = "Test"

## Combined ROC
roc.df = rbind(roc.df,train.roc)
ccs = distinct(roc.df[,3:4])
roc.df$trial = factor(roc.df$trial,levels = c("Train","Test"),
                      label = c(paste0("Train: ",ccs$Comp[ccs$trial=="Train"]),
                                paste0("Test: ",ccs$Comp[ccs$trial!="Train"])) 
)

# ## PLot ROC
# tiff(filename = "Figures/figure1a_relativeAbundance.tiff",width = 2.5,height = 2,units = "in",res = 300)
# ggplot(roc.df,aes(spec,sens,col = trial,group = trial,lty = trial))+
#   geom_line(size = .75)+
#   theme_classic()+
#   xlab("1-Specificity")+
#   ylab("Sensitivity")+
#   ggtitle("RF with Rel. Abundance")+
#   geom_abline(slope = 1,col = "grey")+
#   theme(legend.position = "top",
#         legend.title = element_blank(),plot.title = element_text(hjust = .5,size = 8,face = "bold"),
#         panel.grid = element_blank(),
#         axis.line = element_line(size = .5),
#         legend.text = element_text(size = 6),
#         legend.key.height = unit(0.1, "cm"),
#         legend.margin=margin(-5,0,-15,-10),
#         text = element_text(size = 8),
#         axis.title = element_text(face = "bold",size = 8),
#         axis.text = element_text(size = 8),
#   )+
#   guides(color=guide_legend(nrow=2,byrow=TRUE))
# 
# dev.off()

feat = colnames(df[,-1])
test_ran = data.frame(rDirichlet.acomp(1e4,alpha = 1*rep(1,ncol(df[,-1]))))
colnames(test_ran) = colnames(df[,-1])
testh =m1$models[[1]]
preds = predict.train(testh,test_ran,type = "prob" )
preds.null = predict.train(testh,tbl[,feat],type = "prob" )
test_ran = data.frame(Dataset = "Simulated",Status = "Test",preds ,test_ran)
base = data.frame(Dataset = tbl$Dataset,Status = tbl$Status,preds.null,tbl[,feat])
decison_bound = data.frame(Dataset = xx$Dataset,Status = "decision bound",C1 = 1,C2  = 2,xx[,3:5] )


# ## Plot Ternary
pdf(file = "Figures/fig2_isoporp_RA_ranger.pdf",width = 3 ,height = 3)
ggtern(test_ran,aes(V1,V2,V3))+
  geom_point(aes(col = C1),alpha = .5,pch = 16,size = 1)+
  theme_classic()+
  scale_color_distiller(palette = "RdBu")
dev.off()


test_ran = test_ran %>% 
  rename(Class1_Probability = C1)
base = base %>% 
  rename(Class1_Probability = C1)
decison_bound = decison_bound %>% 
  rename(Class1_Probability = C1)

col_ = if_else(base$Status=="C1","Purple","yellow")
pdf(file = "Figures/fig2_isoprop_RA_ranger.pdf",width = 3 ,height = 3)
ggtern(test_ran,aes(V1,V2,V3,col = Class1_Probability,value = Class1_Probability*1000))+
  geom_point(alpha = .5,pch = 15,size = 3)+
  geom_line(data = decison_bound,aes(x = V1,y = V2,z = V3),size = 1,col  ="black")+
  theme_classic()+
  theme_nolabels()+
  theme_notitles()+
  scale_color_distiller(palette = "RdBu")+
  scale_fill_distiller(palette = "RdBu")+
  geom_point(data = base,aes(V1,V2,V3,shape = Dataset),size  = .75,alpha = .8,col = col_)+
  scale_shape_manual(values = c(4,1))+
  annotate(geom = "text",x = .4,y = .32,z =.18 ,label = "Decision Boundary",angle = 78,size = 2)+
  # annotate(geom = "text",x = .75,y = .15,z =.2 ,label = "Test \nDataset",size = 2,angle = 20)+
  # annotate(geom = "text",x = .08,y = .7,z =.22 ,label = "Train \nDataset",size = 2,angle = 20)+
  theme(legend.position = "none",
        #legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(size = .5),
        plot.margin = margin(-1,-1,-1,-1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        legend.margin=margin(-1,1,-1,-1),
        text = element_text(size = 8),
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
  )  
dev.off()






## Train Model Relative Abundance
## Select Dataset 1
df = tbl %>% 
  filter(Dataset == "Train") %>% 
  select(-Dataset)
df$Status = factor(df$Status)

## Daatset 2
df2 = tbl %>% 
  filter(Dataset == "Test") %>% 
  select(-Dataset)
df2$Status = factor(df2$Status)

## Train RF 
m1 = trainML_Models(trainLRs = df[,-1],testLRs = df2[,-1],
                    ytrain = df[,1],y_test = df2[,1],models = "pls",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

## Train Performance
m1$performance
## Train Performance
classes = as.character(unique(df$Status))
## Train AUC
preds = m1$models[[1]]$pred
train_mroc = pROC::multiclass.roc(preds$obs,preds[,classes])
train.roc = DiCoVarML::rocPlot(train_mroc)
train.roc$trial = "Train"

## Compute AUC
roc_ = pROC::roc(df2$Status,m1$predictionMatrix[,2])
roc_$thresholds
mroc = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,1:2])
roc.df = DiCoVarML::rocPlot(mroc)
roc.df$trial = "Test"

## Combined ROC
roc.df = rbind(roc.df,train.roc)
ccs = distinct(roc.df[,3:4])
roc.df$trial = factor(roc.df$trial,levels = c("Train","Test"),
                      label = c(paste0("Train: ",ccs$Comp[ccs$trial=="Train"]),
                                paste0("Test: ",ccs$Comp[ccs$trial!="Train"])) 
)

# ## PLot ROC
# tiff(filename = "Figures/figure1a_relativeAbundance.tiff",width = 2.5,height = 2,units = "in",res = 300)
# ggplot(roc.df,aes(spec,sens,col = trial,group = trial,lty = trial))+
#   geom_line(size = .75)+
#   theme_classic()+
#   xlab("1-Specificity")+
#   ylab("Sensitivity")+
#   ggtitle("RF with Rel. Abundance")+
#   geom_abline(slope = 1,col = "grey")+
#   theme(legend.position = "top",
#         legend.title = element_blank(),plot.title = element_text(hjust = .5,size = 8,face = "bold"),
#         panel.grid = element_blank(),
#         axis.line = element_line(size = .5),
#         legend.text = element_text(size = 6),
#         legend.key.height = unit(0.1, "cm"),
#         legend.margin=margin(-5,0,-15,-10),
#         text = element_text(size = 8),
#         axis.title = element_text(face = "bold",size = 8),
#         axis.text = element_text(size = 8),
#   )+
#   guides(color=guide_legend(nrow=2,byrow=TRUE))
# 
# dev.off()

feat = colnames(df[,-1])
test_ran = data.frame(rDirichlet.acomp(1e4,alpha = 1*rep(1,ncol(df[,-1]))))
colnames(test_ran) = colnames(df[,-1])
testh =m1$models[[1]]
preds = predict.train(testh,test_ran,type = "prob" )
preds.null = predict.train(testh,tbl[,feat],type = "prob" )
test_ran = data.frame(Dataset = "Simulated",Status = "Test",preds ,test_ran)
base = data.frame(Dataset = tbl$Dataset,Status = tbl$Status,preds.null,tbl[,feat])
decison_bound = data.frame(Dataset = xx$Dataset,Status = "decision bound",C1 = 1,C2  = 2,xx[,3:5] )


# # ## Plot Ternary
# pdf(file = "Figures/fig2_isoporp.pdf",width = 3 ,height = 3)
# ggtern(test_ran,aes(V1,V2,V3))+
#   geom_point(aes(col = C1),alpha = .5,pch = 16,size = 1)+
#   theme_classic()+
#   scale_color_distiller(palette = "RdBu")
# dev.off()


test_ran = test_ran %>% 
  rename(Class1_Probability = C1)
base = base %>% 
  rename(Class1_Probability = C1)
decison_bound = decison_bound %>% 
  rename(Class1_Probability = C1)

col_ = if_else(base$Status=="C1","Purple","yellow")
pdf(file = "Figures/fig2_isoprop_RA_pls.pdf",width = 3 ,height = 3)
ggtern(test_ran,aes(V1,V2,V3,col = Class1_Probability,value = Class1_Probability*1000))+
  geom_point(alpha = .5,pch = 15,size = 3)+
  geom_line(data = decison_bound,aes(x = V1,y = V2,z = V3),size = 1,col  ="black")+
  theme_classic()+
  theme_nolabels()+
  theme_notitles()+
  scale_color_distiller(palette = "RdBu")+
  scale_fill_distiller(palette = "RdBu")+
  geom_point(data = base,aes(V1,V2,V3,shape = Dataset),size  = .75,alpha = .8,col = col_)+
  scale_shape_manual(values = c(4,1))+
  annotate(geom = "text",x = .4,y = .32,z =.18 ,label = "Decision Boundary",angle = 78,size = 2)+
  # annotate(geom = "text",x = .75,y = .15,z =.2 ,label = "Test \nDataset",size = 2,angle = 20)+
  # annotate(geom = "text",x = .08,y = .7,z =.22 ,label = "Train \nDataset",size = 2,angle = 20)+
  theme(legend.position = "none",
        #legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(size = .5),
        plot.margin = margin(-1,-1,-1,-1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        legend.margin=margin(-1,1,-1,-1),
        text = element_text(size = 8),
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
  )  
dev.off()



























# Log ratios --------------------------------------------------------------

### PLR Approach
plr = calcLogRatio(df)
plr2 = calcLogRatio(df2)
m1.plr = trainML_Models(trainLRs = plr[,-1],testLRs = plr2[,-1],
                        ytrain = plr[,1],y_test = plr2[,1],models = "ranger",
                        mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                        numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5)



## Train Performance
m1.plr$performance
classes = as.character(unique(df$Status))
## Train AUC
preds = m1.plr$models[[1]]$pred
train_mroc = pROC::multiclass.roc(preds$obs,preds[,classes])
train.roc = DiCoVarML::rocPlot(train_mroc)
train.roc$trial = "Train"
## Compute AUC
mroc = pROC::multiclass.roc(plr2$Status,m1.plr$predictionMatrix[,1:2])
roc.df = DiCoVarML::rocPlot(mroc)
roc.df$trial = "Test"
## Combined ROC
roc.df = rbind(roc.df,train.roc)
## Combined ROC
roc.df = rbind(roc.df,train.roc)
ccs = unique(roc.df$Comp)
roc.df = rbind(roc.df,train.roc)
ccs = distinct(roc.df[,3:4])
roc.df$trial = factor(roc.df$trial,levels = c("Train","Test"),
                      label = c(paste0("Train: ",ccs$Comp[ccs$trial=="Train"]),
                                paste0("Test: ",ccs$Comp[ccs$trial!="Train"])) 
)

# ## PLot ROC
# tiff(filename = "Figures/figure1a_logRatio.tiff",width = 2.5,height = 2,units = "in",res = 300)
# ggplot(roc.df,aes(spec,sens,col = trial,group = trial,lty = trial))+
#   geom_line(size = .75)+
#   theme_classic()+
#   xlab("1-Specificity")+
#   ylab("Sensitivity")+
#   ggtitle("RF with Log Ratios")+
#   geom_abline(slope = 1,col = "grey")+
#   theme(legend.position = "top",
#         legend.title = element_blank(),plot.title = element_text(hjust = .5,size = 8,face = "bold"),
#         panel.grid = element_blank(),
#         axis.line = element_line(size = .5),
#         legend.text = element_text(size = 6),
#         legend.key.height = unit(0.1, "cm"),
#         legend.margin=margin(-5,0,-15,-10),
#         text = element_text(size = 8),
#         axis.title = element_text(face = "bold",size = 8),
#         axis.text = element_text(size = 8),
#   )+
#   guides(color=guide_legend(nrow=2,byrow=TRUE))
# 
# dev.off()

test_ran = data.frame(rDirichlet.acomp(1e4,alpha = 1*c(1,1,1)))
colnames(test_ran) = colnames(df[,-1])
tt = calcLogRatio(data.frame(Status = "test",test_ran))
testh = m1.plr$models[[1]]
preds = predict.train(testh,tt[,-1],type = "prob" )
preds.null = predict.train(testh,rbind(plr,plr2)[,-1],type = "prob" )
test_ran = data.frame(Status = "Test",preds ,test_ran)
base = data.frame(Status = tbl$Status,preds.null,tbl)
# ## Plot Ternary
ggtern(test_ran,aes(V1,V2,V3))+
  geom_point(aes(col = C1),alpha = .5,pch = 15,size = 2)+
  theme_classic()+
  scale_color_distiller(palette = "RdBu")

decison_bound = decison_bound %>% 
  rename(C1 = Class1_Probability)

col_ = if_else(base$Status=="C1","Purple","yellow")

pdf(file = "Figures/fig2_isoprop_logratio_LR_ranger.pdf",width = 3 ,height = 3)
ggtern(test_ran,aes(V1,V2,V3,col = C1,value = C1*1000))+
  geom_point(alpha = .5,pch = 15,size = 3)+
  geom_line(data = decison_bound,aes(x = V1,y = V2,z = V3),size = 1,col  ="black")+
  theme_classic()+
  theme_nolabels()+
  theme_notitles()+
  scale_color_distiller(palette = "RdBu")+
  scale_fill_distiller(palette = "RdBu")+
  geom_point(data = base,aes(V1,V2,V3,shape = Dataset),size  = 1.25,alpha = 1,col = col_)+
  scale_shape_manual(values = c(4,16))+
  # annotate(geom = "text",x = .4,y = .32,z =.18 ,label = "Decision \nBound",angle = 78,size = 2)+
  # annotate(geom = "text",x = .75,y = .15,z =.2 ,label = "Test \nDataset",size = 2,angle = 20)+
  # annotate(geom = "text",x = .08,y = .7,z =.22 ,label = "Train \nDataset",size = 2,angle = 20)+
  theme(legend.position = "none",
        #legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(size = .5),
        plot.margin = margin(-1,-1,-1,-1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        legend.margin=margin(-1,1,-1,-1),
        text = element_text(size = 8),
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
  )  
dev.off()



































## Train Performance
m1$performance
## Train Performance
classes = as.character(unique(df$Status))
## Train AUC
preds = m1$models[[1]]$pred
train_mroc = pROC::multiclass.roc(preds$obs,preds[,classes])
train.roc = DiCoVarML::rocPlot(train_mroc)
train.roc$trial = "Train"

## Compute AUC
roc_ = pROC::roc(df2$Status,m1$predictionMatrix[,2])
roc_$thresholds
mroc = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,1:2])
roc.df = DiCoVarML::rocPlot(mroc)
roc.df$trial = "Test"

## Combined ROC
roc.df = rbind(roc.df,train.roc)

## PLot ROC
ggplot(roc.df,aes(spec,sens,col = Comp,lty = trial,group = Comp))+
  geom_line(size = 1)+
  theme_bw()+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  geom_abline(slope = 1,col = "grey")+
  theme(legend.position = c(.5,.23),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(size = .5),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.35, "cm"),
        text = element_text(size = 8),
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
  )

feat = colnames(df[,-1])
test_ran = data.frame(rDirichlet.acomp(1e4,alpha = 1*rep(1,ncol(df[,-1]))))
colnames(test_ran) = colnames(df[,-1])
testh =m1$models[[1]]
preds = predict.train(testh,test_ran,type = "prob" )
preds.null = predict.train(testh,df[,feat],type = "prob" )
test_ran = data.frame(ds = "Simulated",Status = "Test",preds ,test_ran)
base = data.frame(ds = "Base",Status = tbl$Status,preds.null,tbl[,feat])
decison_bound = data.frame(ds = "decison",Status = "decision bound",S1 = 1,S2  = 2,xx[,3:5] )


# ## Plot Ternary
ggtern(test_ran,aes(V1,V2,V3))+
  geom_point(aes(col = S1),alpha = .5,pch = 15,size = 2)+
  theme_classic()+
  scale_color_distiller(palette = "RdBu")+
  geom_line(data = xx,aes(x = V1,y = V2,z = V3))

test_ran = test_ran %>% 
  rename(Class1_Probability = S1)
base = base %>% 
  rename(Class1_Probability = S1)
decison_bound = decison_bound %>% 
  rename(Class1_Probability = S1)

col_ = if_else(base$Status=="S1","Purple","yellow")
tiff(filename = "Figures/figure1A.tiff",width = 3,height = 3,units = "in",res = 300)
ggtern(test_ran,aes(V1,V2,V3,col = Class1_Probability,value = Class1_Probability*1000))+
  geom_point(alpha = .5,pch = 15,size = 3)+
  geom_line(data = decison_bound,aes(x = V1,y = V2,z = V3),size = 1,col  ="black")+
  theme_classic()+
  theme_nolabels()+
  scale_color_distiller(palette = "RdBu")+
  scale_fill_distiller(palette = "RdBu")+
  geom_point(data = base,aes(V1,V2,V3),size  = .5,alpha = 1,col = col_)+
  annotate(geom = "text",x = .4,y = .32,z =.18 ,label = "Decision \nBound",angle = 78,size = 2)+
  annotate(geom = "text",x = .75,y = .15,z =.2 ,label = "Test \nDataset",size = 2,angle = 20)+
  annotate(geom = "text",x = .08,y = .7,z =.22 ,label = "Train \nDataset",size = 2,angle = 20)+
  theme(legend.position = c(.85,.8),
        #legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(size = .5),
        plot.margin = margin(-1,-1,-1,-1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(-1,0,-1,-1),
        text = element_text(size = 8),
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
  )  
dev.off()
  


### PLR Approach
plr = calcLogRatio(df)
plr2 = calcLogRatio(df2)
m1.plr = trainML_Models(trainLRs = plr[,-1],testLRs = plr2[,-1],
                        ytrain = plr[,1],y_test = plr2[,1],models = "ranger",
                        mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                        numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5)



## Train Performance
m1.plr$performance
classes = as.character(unique(df$Status))
## Train AUC
preds = m1.plr$models[[1]]$pred
train_mroc = pROC::multiclass.roc(preds$obs,preds[,classes])
train.roc = DiCoVarML::rocPlot(train_mroc)
train.roc$trial = "Train"
## Compute AUC
mroc = pROC::multiclass.roc(plr2$Status,m1.plr$predictionMatrix[,1:2])
roc.df = DiCoVarML::rocPlot(mroc)
roc.df$trial = "Test"
## Combined ROC
roc.df = rbind(roc.df,train.roc)
## PLot ROC
ggplot(roc.df,aes(spec,sens,col = Comp,lty = trial,group = Comp))+
  geom_line(size = 1)+
  theme_bw()+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  geom_abline(slope = 1,col = "grey")+
  theme(legend.position = c(.5,.23),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(size = .5),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.35, "cm"),
        text = element_text(size = 8),
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
  )


test_ran = data.frame(rDirichlet.acomp(1e4,alpha = 1*c(1,1,1)))
colnames(test_ran) = colnames(df[,-1])
tt = calcLogRatio(data.frame(Status = "test",test_ran))
testh = m1.plr$models[[1]]
preds = predict.train(testh,tt[,-1],type = "prob" )
preds.null = predict.train(testh,plr[,-1],type = "prob" )
test_ran = data.frame(Status = "Test",preds ,test_ran)
base = data.frame(Status = tbl$Status,preds.null,tbl)
# ## Plot Ternary
ggtern(test_ran,aes(V1,V2,V3))+
  geom_point(aes(col = S1),alpha = .5,pch = 15,size = 2)+
  theme_classic()+
  scale_color_distiller(palette = "RdBu")+
  geom_line(data = xx,aes(x = V1,y = V2,z = V3))
# 

col_ = if_else(base$Status=="S1","Purple","yellow")
tiff(filename = "Figures/figure1A_logratio.tiff",width = 3,height = 3,units = "in",res = 300)
ggtern(test_ran,aes(V1,V2,V3,col = S1,value = S1*1000))+
  geom_point(alpha = .5,pch = 15,size = 3)+
  geom_line(data = decison_bound,aes(x = V1,y = V2,z = V3),size = 1,col  ="black")+
  theme_classic()+
  theme_nolabels()+
  scale_color_distiller(palette = "RdBu")+
  scale_fill_distiller(palette = "RdBu")+
  geom_point(data = base,aes(V1,V2,V3),size  = .5,alpha = 1,col = col_)+
  annotate(geom = "text",x = .4,y = .32,z =.18 ,label = "Decision \nBound",angle = 78,size = 2)+
  annotate(geom = "text",x = .75,y = .15,z =.2 ,label = "Test \nDataset",size = 2,angle = 20)+
  annotate(geom = "text",x = .08,y = .7,z =.22 ,label = "Train \nDataset",size = 2,angle = 20)+
  theme(legend.position = "none",
        #legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(size = .5),
        plot.margin = margin(-1,-1,-1,-1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(-1,0,-1,-1),
        text = element_text(size = 8),
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
  )  
dev.off()
  
  
  
  
  