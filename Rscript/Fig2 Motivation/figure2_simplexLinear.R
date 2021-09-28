
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
comp =  compositions::rDirichlet.acomp(2,alpha = 200*c(4,2,1));plot(comp);comp = data.frame(comp)[1,]
s1 = data.frame(Dataset = "D1",Status = "C1",comp)

cnames = colnames(comp)
rwnames = rownames(comp)
madison = data.frame()
nsamps = nrow(comp)
alphaRow = which(a==1.2)
lrs.df = data.frame()
a = round(seq(-10,10,by = .1),1)  #<===== Input =======

## Define Decision Boundary
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

x = .5*( (2*madison$V2+madison$V3) / (rowSums(madison[,-2:-1])) )
y = (sqrt(3)/2)*( (madison$V3) / (rowSums(madison[,-2:-1])) )
dec_bound = data.frame(x,y,prob = 0)
dec_bound1 = rbind(dec_bound[-1,],dec_bound[nrow(dec_bound),])
colnames(dec_bound1) = c("x1","y1","prob1")
dec_bound = cbind.data.frame(dec_bound,dec_bound1)


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



tbl = rbind(s2,s3)
tbl$Dataset = if_else(tbl$V1>.33,"Test","Train")








## define grid

# Define Test Grid
ph = data.frame()
incr  =.005
for(a in seq(0,1,by = incr)){
  leftover = 1-a
  b = seq(0,leftover,by =incr)
  c = 1-(b+a)
  ph = rbind(ph,data.frame(a,b,c))
  
}


## convert grid to cartesian
x = .5*( (2*ph$b+ph$c) / (rowSums(ph)) )
y = (sqrt(3)/2)*( (ph$c) / (rowSums(ph)) )
prob = runif(nrow(ph))

cc = data.frame(x,y,prob)
ggplot(cc,aes(x,y,col = prob))+
  geom_point()+
  #geom_point(data = dec_bound,aes(x,y),col = "black")+
  geom_segment(aes(x = x, y = y, xend = x1, yend = y1), data = dec_bound,size = 1,col  ="black")+
  scale_color_distiller(palette = "RdBu")+
  theme_classic()








## Partition Data
df = tbl %>% 
  filter(Dataset == "Train") %>% 
  select(-Dataset)
df$Status = factor(df$Status)

## Daatset 2
df2 = tbl %>% 
  filter(Dataset == "Test") %>% 
  select(-Dataset)
df2$Status = factor(df2$Status)



# Relative Abundance --------------------------------------------------------------

m1 = trainML_Models(trainLRs = df[,-1],testLRs = df2[,-1],
                    ytrain = df[,1],y_test = df2[,1],models = "svmRadial",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

mroc = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,1:2])


test_ran = ph
colnames(test_ran) = colnames(df[,-1])
testh =m1$models[[1]]
preds = predict.train(testh,test_ran,type = "prob" )

cc = data.frame(x,y,prob = preds$C1)
ggplot(cc,aes(x,y,col = prob))+
  geom_point()+
  geom_segment(aes(x = x, y = y, xend = x1, yend = y1), data = dec_bound,size = 1,col  ="black")+
  scale_color_distiller(palette = "RdBu")+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())








# Log Ratios --------------------------------------------------------------


### PLR Approach
plr = calcLogRatio(df)
plr2 = calcLogRatio(df2)
m1.plr = trainML_Models(trainLRs = plr[,-1],testLRs = plr2[,-1],
                        ytrain = plr[,1],y_test = plr2[,1],models = "gbm",
                        mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                        numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5)


## compute ratios
test_ran = ph ## <-------------- Input

colnames(test_ran) = colnames(df[,-1])
tt = calcLogRatio(data.frame(Status = "test",test_ran))
testh = m1.plr$models[[1]]
preds = predict.train(testh,tt[,-1],type = "prob" )
mroc = pROC::multiclass.roc(plr2$Status,m1.plr$predictionMatrix[,1:2])

border = data.frame(start_x = c(0,.5,1),start_y = c(0,sqrt(3)/2,0),
                    end_x = c(.5,1,0),
                    end_yy = c(sqrt(3)/2,0,0),prob = 0)

cc = data.frame(x,y,prob = preds$C1)
ggplot(cc,aes(x,y,col = prob))+
  geom_point()+
  scale_color_distiller(palette = "RdBu")+
  theme_classic()+
  geom_segment(aes(x = x, y = y, xend = x1, yend = y1), data = dec_bound,size = 1,col  ="black")+
  #geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())



