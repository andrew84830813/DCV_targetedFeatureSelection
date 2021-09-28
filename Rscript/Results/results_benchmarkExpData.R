

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
 nm = "cmg_QinN-2014_cirr"
 nm[2] = "cmg_RubelMA-2020_STH"
 nm[3] = "cmg_ZhuF-2020_schizo"
 nm[4] = "cmg_ZellerG_2014_crc"
 
 

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
  ph$col = "lightblue"
  i = which.max(ph$AUC)
  ph$col[i] = "red"
  res.df = rbind(res.df,ph)
}
res = data.frame(res)
results_all$seed_fold = paste0(results_all$Seed,"_",results_all$corrected_fold) 

## Pairwsie comparsion
cx = list( c("DCV-rfRFE","CLR-LASSO"), c("DCV-rfRFE","Coda-LASSO"),
           c("DCV-ridgeEnsemble","CLR-LASSO"), c("DCV-ridgeEnsemble","Coda-LASSO"), 
           c("DCV-ridgeRegression","CLR-LASSO"),c("DCV-ridgeRegression","Coda-LASSO")  )



res1 = results_all %>%
  group_by(Dataset) %>%
  summarise_all(mean) %>% 
  select(Dataset,AUC)
colnames(res1)[2] = "y.position"
dd = results_all %>%
  dplyr::group_by(Dataset) %>%
  rstatix::wilcox_test(data =., AUC ~ Approach,paired = T,comparisons =  cx) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::arrange(dplyr::desc(-p)) %>% 
  filter(p.adj<0.05) %>% 
  left_join(res1)


#tiff(filename =paste0(f_name,".tiff"),width = 4.5,height = 5.5,units = "in",res = 300)
pdf(file = "Figures/benchamrk_ExpData.pdf",width = 7 ,height = 4)
ggplot(results_all,aes(Approach,AUC,shape = Approach))+
  theme_bw()+
  #coord_flip()+
  stat_summary(fun.y = mean, geom = "line",size = 1,col = "black",aes(group =1))+
  #stat_summary(fun.data = mean_se,geom = "errorbar",width = .05)+
  geom_signif(comparisons = cx,tip_length = .05,
              test = "wilcox.test",hjust = -.5,vjust = -1,textsize = 3,
              test.args = list(paired = T),y_position = seq(1,1.25,length.out = 6),
              map_signif_level = function(p) sprintf("p = %.2g", p)
             )+
  geom_line(aes(group =seed_fold),col = "gray")+
  geom_point(col  ="gray",size = 1)+
  scale_shape_manual(values = 21:26)+
  geom_point(data = res.df,aes(Approach,AUC),fill = res.df$col,col = "black",size = 4)+
  facet_wrap(.~Dataset,nrow = 2,
             scales = "free_y"
             )+
  theme(legend.position = "top",
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 8),
        axis.title.y = element_blank(),
        #axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
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

