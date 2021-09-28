toGephi = function(Graph,Name,attribute_metadata=NULL,edge_metadata=NULL){
  el = igraph::as_edgelist(Graph)
  colnames(el) = c("Source","Target")
  el = as.data.frame(el)
  
  edgeWeight = as.matrix(E(Graph)$weight)
  edgeWeight = as.data.frame(edgeWeight)
  edgeWeight = cbind(el,edgeWeight)
  edgeWeight$Ratio = paste0(edgeWeight$Source,"___",edgeWeight$Target)
  colnames(edgeWeight)[3] = "Weight"
  
  if(!is.null(edge_metadata)){
    edgeWeight = left_join(edgeWeight,edge_metadata)
  }
  write_csv(edgeWeight,paste(Name,"_edgeWeights_",".csv",sep= ""))
  
  ddf = (Graph)
  ddf = V(ddf)$name
  ddf = as.data.frame(ddf)
  colnames(ddf) = "KO"
  ex = cbind.data.frame(ID = ddf[,1],Label = ddf[,1])
  if(!is.null(attribute_metadata)){
    ex = left_join(ex,attribute_metadata)
  }
  write_csv(ex,paste(Name,"_attributes_",".csv",sep=""))
  
  
}
