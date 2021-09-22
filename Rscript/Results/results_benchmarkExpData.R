

fnames = dir("Results/")


nm = "cmg_RubelMA-2020_STH"
nm[2] = "cmg_ZhuF-2020_schizo"
nm[3] = "cmg_QinN-2014_cirr"
#nm[4] = "cmg_NielsenHB-2014_ibd"
nm[4] = "cmg_FengQ-2015_crc"
#nm[6] = "cmg_ThomasAM_2019_crc" ## sim
 nm[7] = "cmg_WirbelJ-2018_crc"
# nm[8] = "cmg_YachidaS-2019_crc" ## sim
 nm[9] = "cmg_ZellerG_2014_crc"
# #f_name = "cmg_ZVogtmannE_2016_crc" ## sim
 nm[10] = "cmg_YuJ_2015_crc"

 
 
 
 
 
 nm = "cmg_FengQ-2015_crc"
 nm[2] = "cmg_WirbelJ-2018_crc"
 nm[3] = "cmg_ZellerG_2014_crc"
 nm[4] = "cmg_YuJ_2015_crc"
 
 

results_all = data.frame()
for(f_name in nm){
  bool = str_detect(fnames,paste0(f_name,"_seed"))
  f = fnames[bool]
  results = data.frame()
  for(i in f){
    ph = read_csv(file = paste0("Results/",i))
    results = rbind(results,ph)
  }
  results = separate(results,col = 2,into = c("Dataset","s"),sep = "_seed") %>% 
    dplyr::select(-s)
  ## correct fold mislabeling
  results$corrected_fold = rep(c(rep(1,5),rep(2,5)),5)
  results_all = rbind(results_all,results)
}


res = results_all %>%
  group_by(Approach,Dataset) %>%
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
ggplot(results_all,aes(Approach,AUC))+
  theme_bw()+
  coord_flip()+
  stat_summary(fun.y = mean, geom = "line",size = .75,col = "grey",aes(group =1))+
  geom_point(data = res.df,aes(Approach,AUC),col = res.df$col,size = 3)+
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



