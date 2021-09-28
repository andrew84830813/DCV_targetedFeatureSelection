

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


results$seed_fold = paste0(results$Seed,"_",results$Fold) 
keep_test = c("mbiome_DCV-ridgeRegression")
results = results %>% 
  filter(Approach %in% keep_test)
results1 = results %>% 
  dplyr::select(Approach,seed_fold,Scenario,AUC) %>% 
  spread("Scenario","AUC") %>% 
  mutate(Diff = Empirical-Permuted)
ggplot(results1,aes(Diff))+
  geom_histogram(col = "black",fill = "lightblue")+
  theme_classic()+
  geom_vline(xintercept = mean(results1$Diff),col = "red",lty = "dotdash",size = 1)+
  geom_vline(xintercept = mean(results1$Diff) + (1.96*sd(results1$Diff))/10,col = "red",lty = "dashed",size = 1)+
geom_vline(xintercept = mean(results1$Diff) - (1.96*sd(results1$Diff))/10,col = "red",lty = "dashed",size = 1)
wilcox.test(results1$Empirical,results1$Permuted,paired = T,alternative = "greater")


pdf(file = "Figures/caseStudy_NEC_glmAUC.pdf",width = 2.5 ,height = 3)
ggplot(results,aes(Scenario,AUC,fill = Scenario,label =AUC))+
  theme_bw()+
  #facet_grid(.~Approach)+
  geom_line(aes(group = seed_fold),col = "gray",alpha = .5)+
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
              test.args = list(paired = T) 
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


## compute wilcox test
results1 = results %>% 
  select(seed_fold,Scenario,AUC) %>% 
  spread("Scenario","AUC")
wilcox.test(results1$Empirical,results1$Permuted,paired = T,alternative = "greater")



dd = results %>%
  dplyr::group_by(Approach) %>%
  rstatix::wilcox_test(data =., AUC ~ Scenario,paired = T,alternative = "greater") %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::arrange(dplyr::desc(-p))



