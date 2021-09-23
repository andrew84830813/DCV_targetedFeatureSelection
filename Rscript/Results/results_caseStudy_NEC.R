

fnames = dir("Results/")


## Final Data Sets to include
nm = "NEC_seed"


for(f_name in nm){
  bool = str_detect(fnames,paste0(f_name))
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

n_fun <- function(x){
  return(data.frame(y = mean(x), label = round(mean(x),3) ))
}

#tiff(filename =paste0(f_name,".tiff"),width = 4.5,height = 5.5,units = "in",res = 300)
results$seed_fold
ggplot(results,aes(Scenario,AUC,fill = Scenario,label =AUC))+
  theme_bw()+
  facet_wrap(.~Approach,nrow = 2)+
  geom_line(aes(group = 1),col = "gray",alpha = .5)+
  geom_violin(alpha = .5,)+
  ggsci::scale_color_d3()+
  ggsci::scale_fill_d3()+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.35,
               colour = "black")+
  geom_point(aes(col = Scenario))+
  stat_summary(fun.y = mean, geom = "point",col = "red",size = 2)+
  stat_summary(fun.data = n_fun, geom = "text",size = 4,position = position_nudge(x = .4))+
  stat_summary(fun.data = mean_cl_normal,geom = "errorbar",width = .125)+
  theme(legend.position = "top",
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        axis.title = element_text(size = 8),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        #legend.background = element_rect(colour = "black")
  )





## Pairwsie comparsion
cx = list( c("DCV-rfRFE","CLR-LASSO"), c("DCV-rfRFE","Coda-LASSO"),
           c("DCV-ridgeEnsemble","CLR-LASSO"), c("DCV-ridgeEnsemble","Coda-LASSO"), 
           c("DCV-ridgeRegression","CLR-LASSO"),c("DCV-ridgeRegression","Coda-LASSO")  )


dd = results %>%
  dplyr::group_by(Approach) %>%
  rstatix::wilcox_test(data =., AUC ~ Scenario,paired = T) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::arrange(dplyr::desc(-p))



