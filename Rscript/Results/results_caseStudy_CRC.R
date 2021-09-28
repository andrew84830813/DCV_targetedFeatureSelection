## Top 50
fnames = dir("Results/CRC_top50/")

## Final Data Sets to include
nm = "crcLODO"

for(f_name in nm){
  bool = str_detect(fnames,paste0(f_name,"_seed"))
  f = fnames[bool]
  results = data.frame()
  for(i in f){
    ph = read_csv(file = paste0("Results/CRC_top50/",i))
    results = rbind(results,ph)
  }
}


## Optimal
fnames = dir("Results/")

## Final Data Sets to include
nm = "crcLODO"

for(f_name in nm){
  bool = str_detect(fnames,paste0(f_name,"_seed"))
  f = fnames[bool]
  results = data.frame()
  for(i in f){
    ph = read_csv(file = paste0("Results/",i))
    results = rbind(results,ph)
  }
}


res = results %>%
  group_by(Approach,Scenario,Seed) %>%
  summarise_all(mean)

results = results %>% 
  filter(Approach=="DCV-ridgeRegression")
results$correctedFold = rep(c(1,2,1,2),252)

n_fun <- function(x){
  return(data.frame(y = mean(x), label = round(mean(x),3) ))
}

pdf(file = "Figures/caseStudy_CRC_glmAUC.pdf",width = 2.5 ,height = 3)
ggplot(results,aes(Scenario,AUC,fill = Scenario,label =AUC))+
  theme_bw()+
  #facet_grid(.~Approach)+
  geom_line(aes(group = Seed),col = "gray",alpha = .5)+
  geom_boxplot(alpha = .5,width = .4,outlier.color = NA)+
  scale_color_manual(values = ggsci::pal_lancet()(6)[5:6])+
  scale_fill_manual(values = ggsci::pal_lancet()(6)[5:6])+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.35,
               colour = "black")+
  geom_point(aes(col = Scenario),alpha = .5,shape = 1)+
  geom_signif(comparisons = list(c("Permuted", "Empirical")),tip_length = 0, 
               map_signif_level=F,test = "wilcox.test",
               #test.args = list(paired = T) 
              )+
  stat_summary(fun.y = mean, geom = "point",col = "red",size = 2)+
  stat_summary(fun.data = n_fun, geom = "text",size = 4,position = position_nudge(x = .4))+
  stat_summary(fun.data = mean_cl_normal,geom = "errorbar",width = .125)+
  xlab("Training Label Distribution")+
  theme(legend.position = "none",
        strip.background = element_blank(),
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        axis.title = element_text(size = 8,face = "bold"),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        #legend.background = element_rect(colour = "black")
  )

dev.off()

ggplot(results_all,aes(Approach,number_parts))+
  theme_bw()+
  coord_flip()+
  stat_summary(fun.y = mean, geom = "col",size = 1,col = "black")+
  stat_summary(fun.data = mean_se,geom = "errorbar",width = .5)+
  facet_wrap(.~Dataset,nrow = 2)+
  theme(legend.position = "top",
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 12),
        #axis.title.y = element_blank(),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 12),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        #legend.background = element_rect(colour = "black")
  )






## Pairwsie comparsion
cx = list( c("DCV-rfRFE","CLR-LASSO"), c("DCV-rfRFE","Coda-LASSO"),
           c("DCV-ridgeEnsemble","CLR-LASSO"), c("DCV-ridgeEnsemble","Coda-LASSO"), 
           c("DCV-ridgeRegression","CLR-LASSO"),c("DCV-ridgeRegression","Coda-LASSO")  )


dd = results_all %>%
  dplyr::group_by(Dataset) %>%
  rstatix::wilcox_test(data =., AUC ~ Approach,paired = T,comparisons =  cx) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::arrange(dplyr::desc(-p))



