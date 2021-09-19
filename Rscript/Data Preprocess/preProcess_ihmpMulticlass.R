## Library 
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)


setwd("C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/Archive/dicovarML_project/")
# Read metadata -----------------------------------------------------------
mode.df = data.frame()

## Read data
metadata = data.frame(read_csv(file = "Data/hmp2_metadata (5).csv"))
meta_data = metadata %>% 
  select(External.ID,Participant.ID,week_num,visit_num,is_inflamed,diagnosis,data_type,biopsy_location,SIBDQ.Score,SES.CD.Score,fecalcal,hbi,sccai,fecalcal_ng_ml)
## keep key metadata
metadata = metadata %>% 
  select(External.ID,Participant.ID,diagnosis,data_type)
## data types 
data_type  = unique(metadata$data_type)




# Process Proteomics ------------------------------------------------------

## read data
proteomic = data.frame(read_tsv("Data/HMP2_proteomics_ecs.tsv"))
## remove ungroups in row 1
proteomic = proteomic[-1,]
## gather and spread by features samples
proteomic = gather(proteomic,key = 'External.ID',value = "Counts",2:ncol(proteomic))
proteomic = spread(proteomic,"Gene","Counts")
## merge with metadata
md = metadata %>% 
  filter(External.ID %in% proteomic$External.ID,data_type == "proteomics")
md2 = meta_data%>% 
  filter(External.ID %in% proteomic$External.ID,data_type == "proteomics")
md %>% 
  group_by(Participant.ID,diagnosis) %>% 
  summarise(n = n()) %>% 
  group_by(diagnosis) %>% 
  summarise(n  = n())
proteomic = right_join(md,proteomic)
md = proteomic[,1:3]
df = proteomic[,c(-1,-2,-4)]
colnames(df)[1] = "Status"
# df$Status = if_else(df$Status!="nonIBD","IBD","nonIBD")
# md$diagnosis = if_else(md$diagnosis!="nonIBD","IBD","nonIBD")
rownames(df) = md$External.ID
## remove samples wit all 0's
bool = !rowSums(df[,-1])==0
md = md[bool,1:3]
df = df[bool,]
## Store data
obj = list(data = df,md = md,all_md = md2)
save(obj,file = "C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/DiCoVarFS_project/Output/ihmp/proteomics_data.Rda")
mode.df = rbind(mode.df,data.frame(External.ID = md$External.ID,data_type = "proteomics" ))




# Metabolites -------------------------------------------------------------

## read data
proteomic = data.frame(read_csv("Data/iHMP_metabolomics.csv"))
## remove ungroups in row 1
proteomic = proteomic %>% 
  filter(!is.na(Metabolite))

modes = unique(proteomic$Method)

for(i in 1:length(modes)){
  ph = proteomic %>% 
    filter(Method==modes[i]) 
  ph = ph[,c(-5:-1,-7)]
  ph = ph %>% 
    group_by(Metabolite) %>% 
    summarise_all(.funs = sum)
  ## gather and spread by features samples
  ph = gather(ph,key = 'External.ID',value = "Counts",2:ncol(ph))
  ph = spread(ph,"Metabolite","Counts",fill = 0)
  ph = data.frame(ph)
  sum(is.na(rowSums(ph[,-1])))
  ## merge with metadata
  md = metadata %>% 
    filter(External.ID %in% ph$External.ID,data_type == "metabolomics")
  md2 = meta_data %>% 
    filter(External.ID %in% ph$External.ID,data_type == "metabolomics")
  md %>% 
    group_by(Participant.ID,diagnosis) %>% 
    summarise(n = n()) %>% 
    group_by(diagnosis) %>% 
    summarise(n  = n())
  ph = right_join(md,ph)
  ## remove samples wit all 0's
  bool = !rowSums(ph[,c(-1,-2,-3,-4)])==0
  md = md[bool,1:3]
  df = ph[bool,c(-1,-2,-4)]
  ## Final Prep
  colnames(df)[1] = "Status"
  # df$Status = if_else(df$Status!="nonIBD","IBD","nonIBD")
  # md$diagnosis = if_else(md$diagnosis!="nonIBD","IBD","nonIBD")
  rownames(df) = md$External.ID
  ## Store data
  obj = list(data = df,md = md,all_md = md2)
  fname = paste0("C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/DiCoVarFS_project/Output/ihmp/metabolomics_",modes[i],"_data.Rda")
  save(obj,file = fname)
  mode.df = rbind(mode.df,data.frame(External.ID = md$External.ID,data_type = modes[i]))
  
}


# Meta Transc  -------------------------------------------------------------

## read data
proteomic = data.frame(read_tsv("Data/hmp2_mtx_pathabundance.tsv"))
proteomic = separate(proteomic,col = 1,into = c("Path","Species"),sep = "\\|")
proteomic = proteomic[-2:-1,]
bool  =is.na(proteomic$Species)
proteomic = proteomic[bool,]

## remove ungroups in row 1
proteomic = proteomic[,-2]
proteomic = proteomic %>% 
  group_by(Path) %>% 
  summarise_all(.funs = sum)
## gather and spread by features samples
proteomic = gather(proteomic,key = 'External.ID',value = "Counts",2:ncol(proteomic))
proteomic = spread(proteomic,"Path","Counts")
id = str_split(proteomic$External.ID,pattern = "_pathabundance_cpm",n = 2,simplify = T)[,1]
proteomic$External.ID = id
## merge with metadata
md = metadata %>% 
  filter(External.ID %in% proteomic$External.ID,data_type == "metatranscriptomics")
md2 = meta_data %>% 
  filter(External.ID %in% proteomic$External.ID,data_type == "metatranscriptomics")
md %>% 
  group_by(Participant.ID,diagnosis) %>% 
  summarise(n = n()) %>% 
  group_by(diagnosis) %>% 
  summarise(n  = n())
proteomic = right_join(md,proteomic)
md = proteomic[,1:3]
df = proteomic[,c(-1,-2,-4)]
colnames(df)[1] = "Status"
# df$Status = if_else(df$Status!="nonIBD","IBD","nonIBD")
# md$diagnosis = if_else(md$diagnosis!="nonIBD","IBD","nonIBD")
rownames(df) = md$External.ID
## remove samples wit all 0's
bool = !rowSums(df[,-1])==0
md = md[bool,1:3]
df = df[bool,]
## remove samples wit all 0's
bool = !is.na(df$Status)
md = md[bool,1:3]
df = df[bool,]
## Store data
obj = list(data = df,md = md,all_md = md2)
save(obj,file = "C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/DiCoVarFS_project/Output/ihmp/metatranscriptomics_data.Rda")
mode.df = rbind(mode.df,data.frame(External.ID = md$External.ID,data_type = "metatranscriptomics" ))





#  mgx Abundance ----------------------------------------------------------

## read data
proteomic = data.frame(read_tsv("Data/hmp2_mgx_pathabundance.tsv"))
proteomic = separate(proteomic,col = 1,into = c("Path","Species"),sep = "\\|")
proteomic = proteomic[-2:-1,]
bool  =is.na(proteomic$Species)
proteomic = proteomic[bool,]

## remove ungroups in row 1
proteomic = proteomic[,-2]
proteomic = proteomic %>% 
  group_by(Path) %>% 
  summarise_all(.funs = sum)
## gather and spread by features samples
proteomic = gather(proteomic,key = 'External.ID',value = "Counts",2:ncol(proteomic))
proteomic = spread(proteomic,"Path","Counts")
id = str_split(proteomic$External.ID,pattern = "_pathabundance_cpm",n = 2,simplify = T)[,1]
proteomic$External.ID = id
## merge with metadata
md = metadata %>% 
  filter(External.ID %in% proteomic$External.ID,data_type == "metagenomics")
md2 = meta_data %>% 
  filter(External.ID %in% proteomic$External.ID,data_type == "metagenomics")
md %>% 
  group_by(Participant.ID,diagnosis) %>% 
  summarise(n = n()) %>% 
  group_by(diagnosis) %>% 
  summarise(n  = n())
proteomic = right_join(md,proteomic)
md = proteomic[,1:3]
df = proteomic[,c(-1,-2,-4)]
colnames(df)[1] = "Status"
#df$Status = if_else(df$Status!="nonIBD","IBD","nonIBD")
#md$diagnosis = if_else(md$diagnosis!="nonIBD","IBD","nonIBD")
rownames(df) = md$External.ID
## remove samples wit all 0's
bool = !rowSums(df[,-1])==0
md = md[bool,1:3]
df = df[bool,]
## remove samples wit all 0's
bool = !is.na(df$Status)
md = md[bool,1:3]
df = df[bool,]
## Store data
obj = list(data = df,md = md,all_md = md2)
save(obj,file = "C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/DiCoVarFS_project/Output/ihmp/mgxFunctions_data.Rda")
mode.df = rbind(mode.df,data.frame(External.ID = md$External.ID,data_type = "mgx_functions" ))





# 
# 
# #  mgx tax ----------------------------------------------------------
# 
## read data
proteomic = data.frame(read_tsv("Data/hmp2_mgx_taxonomy.tsv"))
proteomic = separate(proteomic,col = 1,into = c("K","P","C","O","F","Genus","Species"),sep = "\\|")
## remove ungroups in row 1
proteomic = proteomic[,-6:-1]
proteomic = proteomic %>%
  filter(!is.na(Species)) %>% 
  group_by(Species) %>%
  summarise_all(.funs = sum) 
## gather and spread by features samples
proteomic = gather(proteomic,key = 'External.ID',value = "Counts",2:ncol(proteomic))
proteomic = spread(proteomic,"Species","Counts")
id = str_split(proteomic$External.ID,pattern = "_profile",n = 2,simplify = T)[,1]
proteomic$External.ID = id
colnames(proteomic) = str_replace_all(colnames(proteomic),pattern = "\\[",replacement = "")
colnames(proteomic) = str_replace_all(colnames(proteomic),pattern = "\\]",replacement = "")
## merge with metadata
md = metadata %>% 
  filter(External.ID %in% proteomic$External.ID,data_type == "metagenomics")
md2 = meta_data %>% 
  filter(External.ID %in% proteomic$External.ID,data_type == "metagenomics")
md %>% 
  group_by(Participant.ID,diagnosis) %>% 
  summarise(n = n()) %>% 
  group_by(diagnosis) %>% 
  summarise(n  = n())
proteomic = right_join(md,proteomic)
md = proteomic[,1:3]
df = proteomic[,c(-1,-2,-4)]
colnames(df)[1] = "Status"
# df$Status = if_else(df$Status!="nonIBD","IBD","nonIBD")
# md$diagnosis = if_else(md$diagnosis!="nonIBD","IBD","nonIBD")
rownames(df) = md$External.ID
## remove samples wit all 0's
bool = !rowSums(df[,-1])==0
md = md[bool,1:3]
df = df[bool,]
## remove samples wit all 0's
bool = !is.na(df$Status)
md = md[bool,1:3]
df = df[bool,]
## Store data
obj = list(data = df,md = md,all_md = md2)
save(obj,file = "C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/DiCoVarFS_project/Output/ihmp/mgxTaxonomy_data.Rda")
mode.df = rbind(mode.df,data.frame(External.ID = md$External.ID,data_type = "mgx_taxonomy" ))







### mulit omic integration
metadata = data.frame(read_csv(file = "Data/hmp2_metadata (5).csv"))
m = metadata %>% 
  select(External.ID,Participant.ID,week_num,visit_num,is_inflamed,
         diagnosis,data_type,biopsy_location,SIBDQ.Score,SES.CD.Score,fecalcal,hbi,sccai,fecalcal_ng_ml)
m = m %>% 
  group_by(External.ID,Participant.ID,diagnosis) %>% 
  summarise_all(.funs = mean)

mode.df1 = mode.df 
mode.df1$ct = 1
mode.df1 = spread(mode.df1,"data_type","ct",fill = 0)
mode.df1$modes = rowSums(mode.df1[,-1])
table(mode.df1$modes)

mm = left_join(mode.df1,m) %>% 
  filter(modes==8)
table(mm$diagnosis)

n_distinct(mm$Participant.ID)
mm %>% 
  group_by(Participant.ID,diagnosis) %>% 
  summarise(n = n()) %>% 
  group_by(diagnosis) %>% 
  summarise(n = n())
keep_samples = mm
save(keep_samples,file = "C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/DiCoVarFS_project/Output/ihmp/multiomic_keep.Rda")




