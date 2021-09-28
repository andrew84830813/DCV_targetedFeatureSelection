
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


## Define Decision Boundary
LR = mean(log(tbl$V1/tbl$V3))
xx = data.frame(Dataset = "",Status = "Expected Decison Bound",a = seq(.01,1,length.out = 1000) )
xx$c = xx$a/exp(LR)
xx = xx %>% 
  dplyr::mutate(S = a+c) %>% 
  filter(S<=1) %>% 
  dplyr::select(-S) %>% 
  dplyr::mutate(b = 1 - (a + c)) %>% 
  dplyr::select(a,b,c) 
x = .5*( (2*xx$b+xx$c) / (rowSums(xx)) )
y = (sqrt(3)/2)*( (xx$c) / (rowSums(xx)) )
dec_bound = data.frame(x,y,prob = 0)



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
  geom_point(data = dec_bound,aes(x,y),col = "black")+
  scale_color_distiller(palette = "RdBu")+
  theme_classic()


# cc = data.frame(x = round(x*1000),y = round(y*1000))
# cc = distinct(cc)
# cc$prob = runif(nrow(cc))
# xx = dec_bound*1000
# ggplot(cc,aes(x,y,fill = prob))+
#   geom_tile()+
#   scale_fill_distiller(palette = "RdBu")+
#   geom_line(data = xx,aes(x,y),col = "red")+
#   theme_classic()




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
                    ytrain = df[,1],y_test = df2[,1],models = "pls",
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
  geom_line(data = dec_bound,aes(x,y),col = "black")+
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
                        ytrain = plr[,1],y_test = plr2[,1],models = "glmnet",
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
  geom_line(data = dec_bound,aes(x,y),col = "black")+
  scale_color_distiller(palette = "RdBu")+
  theme_classic()+
  geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())

