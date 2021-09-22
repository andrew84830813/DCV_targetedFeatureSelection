
## data locations
fnames = (dir("Output/CRC_exp/"))


crcCombined = data.frame()
for(f in fnames){
  load(paste0("Output/CRC_exp/",f))
  df = expResults$taxa.df
  df = data.frame(sample_id = rownames(df),df)
  rdf = tidyr::gather(df,"ID","Count",3:ncol(df))
  
  cname.table = expResults$taxa_info %>% 
    dplyr::rename(ID = species)
  cname.table$ID = stringr::str_replace(cname.table$ID,pattern = " ",replacement = "\\.")
  rdf = left_join(rdf,cname.table)
  rdf = subset(rdf,select = c("sample_id","Status","ID","Count"))
  rdf = rdf %>% 
    dplyr::group_by(sample_id,ID,Status) %>% 
    dplyr::summarise(count = sum(Count) ) 
  fff = stringr::str_split(f,pattern = "curatedMetaGenome_",simplify = T)
  fff = stringr::str_split(fff[,2],pattern = "\\.",simplify = T)[,1]
  rdf = data.frame(dataset = fff,rdf)
  crcCombined = rbind(crcCombined,rdf)
}

df = tidyr::spread(crcCombined,"ID","count",fill = 0)
md = df[,1:3]
df = data.frame(df[,-2:-1],row.names = df$sample_id)
bool = df$Status %in% c("CRC","control")
df = df[bool,]
md = md[bool,]
table(df$Status)
table(md$dataset)/sum(table(md$dataset))

## total counts
rs = rowSums(df[,-1])
hist(rs)
summary(rs)
unique(md$dataset)

obj = list(data = df,md = md)
save(obj,file = "Output/LODO/crc_data.Rda")


