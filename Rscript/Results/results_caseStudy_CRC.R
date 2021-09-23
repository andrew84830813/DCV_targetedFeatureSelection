

fnames = dir("Results/")


# nm = "cmg_RubelMA-2020_STH"
# nm[2] = "cmg_ZhuF-2020_schizo"
# nm[3] = "cmg_QinN-2014_cirr"
# #nm[4] = "cmg_NielsenHB-2014_ibd"
# nm[4] = "cmg_FengQ-2015_crc"
# #nm[6] = "cmg_ThomasAM_2019_crc" ## sim
#  nm[7] = "cmg_WirbelJ-2018_crc"
# # nm[8] = "cmg_YachidaS-2019_crc" ## sim
#  nm[9] = "cmg_ZellerG_2014_crc"
# # #f_name = "cmg_ZVogtmannE_2016_crc" ## sim
#  nm[10] = "cmg_YuJ_2015_crc"
# 
#  



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

n_fun <- function(x){
  return(data.frame(y = mean(x), label = round(mean(x),3) ))
}

#tiff(filename =paste0(f_name,".tiff"),width = 4.5,height = 5.5,units = "in",res = 300)
ggplot(res,aes(Scenario,AUC,fill = Scenario,label =AUC))+
  theme_bw()+
  facet_grid(.~Approach)+
  geom_line(aes(group = Seed),col = "gray",alpha = .5)+
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



