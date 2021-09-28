
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


# Load Helper Functions  ---------------------------------------------------------------
fnames = dir("Helper_Functions/")
for(f in fnames){
  source(paste0("Helper_Functions/",f))
}


clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)





# Sim Dataset -------------------------------------------------------------

## Shift Data Within Simplex
## D1
comp =  compositions::rDirichlet.acomp(2,alpha = 200*c(4,2,1));plot(comp);comp = data.frame(comp)
s1 = data.frame(Dataset = "D1",Status = "C1",comp)
# comp =  compositions::rDirichlet.acomp(100,alpha = 20*c(3,1,1));plot(comp);comp = data.frame(comp)
# s2 = data.frame(Dataset = "D1",Status = "C2",comp)
# d1 = rbind(s1,s2)


cnames = colnames(comp)
rwnames = rownames(comp)
madison = data.frame()
nsamps = nrow(comp)
alphaRow = which(a==1.2)
lrs.df = data.frame()
a = round(seq(-10,10,by = .1),1)  #<===== Input =======

for(j in 1:nsamps){
  x1 = comp[j,]
  c.df = data.frame()
  label = c()
  i = 1
  for(A in a){
    ph = data.frame((clo(x1^A)))
    c.df = rbind(c.df,ph)
    label[i] = "black"
    i = i+1
  }
  # lrs = calcLogRatio(data.frame(Status = paste0("p",j),c.df))
  # lrs.df = rbind(lrs.df,lrs)#

  #madison = rbind(madison,data.frame(c.df[alphaRow,]))
  madison = rbind(madison,data.frame(Dataset = j,Status = a,c.df))
  #ADD Indv Points
  colnames(c.df) = cnames
  rownames(c.df) = paste("A",as.character(round(a,digits = 3)),sep = "_") #with alpha values
  labels = str_split(rownames(c.df),pattern = "_",simplify = T,n = 2)
  labels = as.factor(labels[,1])
  #plot
  ac = acomp(c.df)
  if(j==1){
    plot(ac,col = labels,"l")
    plot(acomp(comp),add = T,pch = 18,col = "red")
    #plot(compMean,col = "steelblue",add = T,pch = 8,cex = 2)
    #plot(acomp(baryCenter),col = "green",add = T,cex= 2,pch = 16)
  }else{
    plot(ac,col = labels,"l",add = T)
  }
}

ggtern(s1,aes(x = V1,y = V2,z = V3))+
  geom_point(aes(col = Status),size = 1,alpha = .5)+
  geom_line(data = madison,aes(x = V1,y = V2,z = V3,group = Dataset))+
  theme_classic()

alpha = runif(n = 10,min = 0,max = .125)
c = simplexDataAugmentation::linearShift(data.frame(madison[,-2:-1]),
                                         a = alpha,alpha_n = alpha,pertubationDirection = clo(c(100,1,1)),
                                         directionType = "centroid",evComp = 1)
c = c$linAdjustedData
s2 = data.frame(Dataset = "d3",Status = "C1",sample_n(c,size = 200,replace = F))
ggtern(s1,aes(x = V1,y = V2,z = V3))+
  #geom_point(aes(col = Status),size = 1,alpha = .5)+
  geom_point(data = s2,aes(col = Status),size = 1,alpha = .5)+
  geom_line(data = madison,aes(x = V1,y = V2,z = V3,group = Dataset))+
  theme_classic()


alpha = .1
c = simplexDataAugmentation::linearShift(data.frame(s2[,-2:-1]),
                                         a = alpha,alpha_n = alpha,pertubationDirection = clo(c(1,100,1)),
                                         directionType = "centroid",evComp = 1)
c = c$linAdjustedData
s3 = data.frame(Dataset = "d3",Status = "C2",c)
ggtern(s1,aes(x = V1,y = V2,z = V3))+
  #geom_point(aes(col = Status),size = 1,alpha = .5)+
  geom_point(data = s2,aes(col = Status),size = 1,alpha = .5)+
  geom_point(data = s3,aes(col = Status),size = 1,alpha = .5,col = "blue")+
  geom_line(data = madison,aes(x = V1,y = V2,z = V3,group = Dataset))+
  theme_classic()

tbl = rbind(s2,s3)
tbl$Dataset = if_else(tbl$V1>.33,"Test","Train")


compMean =as.numeric( mean.acomp(acomp(tbl[,-2:-1])))
compMean = data.frame(Dataset = "mean",Status = "mean",V1 = compMean[1],V2 =compMean[2],V3=compMean[3])
t = rbind(tbl,compMean)
c = simplexDataAugmentation::linearShift(data.frame(t[,-2:-1]),
                                         a = seq(-10,10,.5),alpha_n = 1,
                                         directionType = "eigen",evComp = 1)
c = c$allShiftData
c =  data.frame(Dataset = "mean",Status = "mean",c[c$sampleNum==401,])


## Overall Data
ggtern(tbl[tbl$Dataset=="Train",],aes(x = V1,y = V2,z = V3,col = Status))+
  geom_point(size = 1,alpha = .5)+
  theme_classic()
ggtern(tbl[tbl$Dataset=="Test",],aes(x = V1,y = V2,z = V3,col = Status))+
  geom_point(size = 1,alpha = .5)+
  theme_classic()
ggtern(tbl,aes(x = V1,y = V2,z = V3,col = Status, shape = Dataset))+
  geom_point(size = 1,alpha = .5)+
  theme_classic()+
  geom_line(data =c, aes(x = V1,y = V2,z = V3))
  























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
                    ytrain = df[,1],y_test = df2[,1],models = "glmnet",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 10
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
# tiff(filename = "Figures/figure1bc_relativeAbundanceROC.tiff",width = 2.5,height = 2,units = "in",res = 300)
# ggplot(roc.df,aes(spec,sens,col = trial,group = trial,lty = trial))+
#   geom_line(size = .75)+
#   theme_classic()+
#   xlab("1-Specificity")+
#   ylab("Sensitivity")+
#   ggtitle("PLS with Rel. Abundance")+
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
# dev.off()

feat = colnames(df[,-1])
test_ran = data.frame(rDirichlet.acomp(1e4,alpha = 1*rep(1,ncol(df[,-1]))))
colnames(test_ran) = colnames(df[,-1])
testh =m1$models[[1]]
preds = predict.train(testh,test_ran,type = "prob" )
preds.null = predict.train(testh,tbl[,feat],type = "prob" )
test_ran = data.frame(Dataset = "Simulated",Status = "Test",preds ,test_ran)
base = data.frame(Dataset = tbl$Dataset,Status = tbl$Status,preds.null,tbl[,feat])
decison_bound = data.frame(Dataset = "decison",Status = "decision bound",C1 = 1,C2  = 2,c )


# ## Plot Ternary
ggtern(test_ran,aes(V1,V2,V3))+
  geom_point(aes(col = C1),alpha = .5,pch = 15,size = 2)+
  theme_classic()+
  scale_color_distiller(palette = "RdBu")

test_ran = test_ran %>% 
  rename(Class1_Probability = C1)
base = base %>% 
  rename(Class1_Probability = C1)
decison_bound = decison_bound %>% 
  rename(Class1_Probability = C1)
col_ = if_else(base$Status=="C1","Purple","yellow")

pdf(file = "Figures/fig2_simplexLinear_RA_ranger.pdf",width = 3 ,height = 3)
ggtern(test_ran,aes(V1,V2,V3,col = Class1_Probability,value = Class1_Probability*1000))+
  geom_point(alpha = .5,pch = 15,size = 3)+
  geom_line(data = decison_bound,aes(x = V1,y = V2,z = V3),size = 1,col  ="black")+
  theme_classic()+
  theme_nolabels()+
  theme_notitles()+
  scale_color_distiller(palette = "RdBu")+
  scale_fill_distiller(palette = "RdBu")+
  #geom_point(data = base,aes(V1,V2,V3,shape = Dataset),size  = .75,alpha = .8,col = col_)+
  scale_shape_manual(values = c(4,1))+
  annotate(geom = "text",x = .4,y = .32,z =.18 ,label = "Decision Boundary",angle = 78,size = 2)+
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



### PLR Approach
plr = calcLogRatio(df)
plr2 = calcLogRatio(df2)
m1.plr = trainML_Models(trainLRs = plr[,-1],testLRs = plr2[,-1],
                        ytrain = plr[,1],y_test = plr2[,1],models = "glmnet",
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

## PLot ROC
tiff(filename = "Figures/figure1b_logRatioROC.tiff",width = 2.5,height = 2,units = "in",res = 300)
ggplot(roc.df,aes(spec,sens,col = trial,group = trial,lty = trial))+
  geom_line(size = .75)+
  theme_classic()+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  ggtitle("PLS with Log Ratios")+
  geom_abline(slope = 1,col = "grey")+
  theme(legend.position = "top",
        legend.title = element_blank(),plot.title = element_text(hjust = .5,size = 8,face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(size = .5),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.1, "cm"),
        legend.margin=margin(-5,0,-15,-10),
        text = element_text(size = 8),
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
  )+
  guides(color=guide_legend(nrow=2,byrow=TRUE))

dev.off()

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
tiff(filename = "Figures/figure1b_logratioTern.tiff",width = 3,height = 3,units = "in",res = 300)
ggtern(test_ran,aes(V1,V2,V3,col = C1,value = C1*1000))+
  geom_point(alpha = .5,pch = 15,size = 3)+
  geom_line(data = decison_bound,aes(x = V1,y = V2,z = V3),size = 1,col  ="black")+
  theme_classic()+
  theme_nolabels()+
  theme_notitles()+
  scale_color_distiller(palette = "RdBu")+
  scale_fill_distiller(palette = "RdBu")+
  #geom_point(data = base,aes(V1,V2,V3,shape = Dataset),size  = 1.25,alpha = 1,col = col_)+
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




