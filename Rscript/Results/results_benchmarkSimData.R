library(dplyr)
library(ggplot2)

fnames = dir("Results/")




## Final Data Sets to include
nm = "16S_meanShift"
nm[2] = "WGS_meanShift"


results_all = data.frame()
for(f_name in nm){
  bool = str_detect(fnames,f_name)
  f = fnames[bool]
  results = data.frame()
  for(i in f){
    ph = read_csv(file = paste0("Results/",i))
    ph$corrected_fold = rep(c(rep(1,5),rep(2,5)),5)
    results = rbind(results,ph)
  }
  results_all= rbind(results_all,results)
}
results_all = tidyr::separate(results_all,col = 2,into = c("Dataset","shift","Permute","S"),sep = "_")
  



res = results_all %>%
  group_by(Scenario,Dataset,shift_parm) %>%
  summarise_all(mean)
ds = unique(res$Dataset)
res.df = data.frame()
for(d in ds){
  ph = res %>% 
    filter(Dataset==d)
  ph$col = "black"
  i = which.max(ph$AUC)
  ph$col[i] = "red"
  res.df = rbind(res.df,ph)
}
res = data.frame(res)

#tiff(filename =paste0(f_name,".tiff"),width = 4.5,height = 5.5,units = "in",res = 300)

results_all1 = results_all %>% 
  filter(Scenario=="Empirical")
results_all2 = results_all %>% 
  filter(Scenario!="Empirical")
ggplot(results_all1,aes(shift_parm,AUC,col = Approach,shape = Approach))+
  theme_bw()+
  stat_summary(fun.y = mean, geom = "line",size = .75,aes(group =Approach))+
  stat_summary(fun.y = mean, geom = "point",size = 3,aes(group =Approach))+
  stat_summary(fun.data = mean_cl_normal,geom = "errorbar",width = .25)+
  facet_wrap(.~Dataset,nrow = 1)+
  theme(legend.position = "top",panel.grid = element_blank(),
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


f <- function(x) {
  r <- quantile(x, probs = c(0.10, 0.25, 0.5, 0.75, 0.90))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

sets = data.frame(shift_parm = 1:6,shift_per = (seq(1,1.3,length.out = 6)-1)*100)
results_all = left_join(results_all,sets)


pdf(file = "Figures/benchamrk_SimData.pdf",width = 7 ,height = 4)
ggplot(results_all,aes(shift_per,AUC,col = Scenario,shape =Scenario))+
  theme_bw()+
  #geom_point()+
  ggsci::scale_color_d3()+
  xlab("% Mean Shift")+
  #scale_x_continuous(breaks = unique(results_all$shift_parm))+
  stat_summary(fun.y = mean, geom = "line",size = .75)+
  stat_summary(fun.y = mean, geom = "point",size = 2)+
  stat_summary(fun.data = ggplot2::mean_cl_normal,geom = "errorbar",width = .25)+
  facet_grid(Dataset~Approach)+
  theme(legend.position = "top",panel.grid = element_blank(),
        strip.background = element_blank(),strip.text = element_text(face = "bold",size = 8),
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
dev.off()







## Pairwsie comparsion
cx = list( c("DCV-rfRFE","CLR-LASSO"), c("DCV-rfRFE","Coda-LASSO"),
           c("DCV-ridgeEnsemble","CLR-LASSO"), c("DCV-ridgeEnsemble","Coda-LASSO"), 
           c("DCV-ridgeRegression","CLR-LASSO"),c("DCV-ridgeRegression","Coda-LASSO")  )


dd = results_all1 %>%
  dplyr::group_by(Scenario,Dataset,shift_parm) %>%
  rstatix::wilcox_test(data =., AUC ~ Approach,paired = T,comparisons =  cx) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::arrange(dplyr::desc(-p))



