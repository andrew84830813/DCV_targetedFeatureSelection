

## Scenario-1 = Sim from 16S Data with mean shift
g1 = 100
g2 = 100
seed_ = 1
sparsePercent = .9



mdwgs = readRDS("Output/16sModel.RDS");
set.seed(seed_);
dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts));
dat = sample_n(dat,size = g1+g2,replace = F);
labels = sample(c(rep("S1",g1),rep("S2",g2)));
dat = data.frame(Status = labels,dat);
fname = "16S_meanShift";
#process shift
procData = processCompData(dat,minPrevalence = sparsePercent);
dat = procData$processedData;
impFact = procData$impFactor;
y = dat[,-1];
bool = colSums(y)==0;
y = y[,!bool];
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact));

dat = simFromExpData.largeMeanShft(raMatrix = dat[,-1],
                                   n1 = g1,n2 = g2,
                                   featureShiftPercent =  1.3,
                                   perFixedFeatures = .9)





mdwgs = readRDS("Output/wgsModel.RDS");
set.seed(seed_);
dat = data.frame(t(zinbwave::zinbSim(mdwgs)$counts));
dat = sample_n(dat,size = g1+g2,replace = F);
labels = sample(c(rep("S1",g1),rep("S2",g2)));
dat = data.frame(Status = labels,dat);
fname = "WGS_meanShift";
#process shift
procData = processCompData(dat,minPrevalence = sparsePercent);
dat = procData$processedData;
impFact = procData$impFactor;
y = dat[,-1];
bool = colSums(y)==0;
y = y[,!bool];
dat = data.frame(Status = dat[,1],fastImputeZeroes( compositions::clo(y),impFactor = impFact));

dat = simFromExpData.largeMeanShft(raMatrix = dat[,-1],
                                   n1 = g1,n2 = g2,
                                   featureShiftPercent =  1.15,
                                   perFixedFeatures = .9)




df = dat
lrs = calcLogRatio(dat)
rt = prcomp(lrs[,-1])
coords.df = data.frame(Group = lrs[,1],rt$x)
ggplot(coords.df,aes(PC1,PC2,col = Group))+
  geom_point(aes(shape = Group,fill  = Group),alpha = .6,size = 3)+
  theme_bw()+
  ggtitle("PCA of all pairwise logratios")+
  theme(legend.position = "top",
        axis.title = element_text(face = "bold",size = 20),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),
        legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))
