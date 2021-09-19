

getCuratedMetagenDataV3_0 = function(fname){
  ds = paste0(fname,".relative_abundance")
  ww = curatedMetagenomicData::curatedMetagenomicData(ds, dryrun = FALSE,counts = T)
  ## get exp obj
  ww = ww[[1]]
  ##get taxa info
  ttt = SummarizedExperiment::rowData(ww)
  ## get count matrix
  count.df = SummarizedExperiment::assay(ww)
  count.df = data.frame(Taxa = ttt$species,count.df)
  count.df = tidyr::gather(count.df,"sampleID","Count",2:ncol(count.df))
  count.df = tidyr::spread(count.df,key = "Taxa",value = "Count",fill = 0)
  ## get patient metadata
  pt.metadata =  as.data.frame(colData(ww))
  pt.metadata = data.frame(sampleID = rownames(pt.metadata),pt.metadata)
  ## merge study condition
  md = subset(pt.metadata,select = c("sampleID","study_condition"))
  md$sampleID = stringr::str_replace_all(string = md$sampleID,pattern = "-",replacement = "\\.")
  df = dplyr::left_join(md,count.df)
  df = data.frame(Status = df$study_condition,df[,-2:-1],row.names = df$sampleID)
  expResults = list(taxa.df = df,metadata = pt.metadata,taxa_info = data.frame(ttt))
  expResults
}



