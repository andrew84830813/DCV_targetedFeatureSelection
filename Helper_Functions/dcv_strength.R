dcv_strength = function (dcvScores) 
{
  Str = NULL
  dcvScores$rowmean[dcvScores$rowmean < 0] = 0
  trainData.md = caret::preProcess(data.frame(dcvScores$rowmean), 
                                   method = "range", rangeBounds = c(0, 1))
  scaledScore = stats::predict(trainData.md, data.frame(dcvScores$rowmean))
  el = data.frame(Ratio = dcvScores$Ratio, Score = scaledScore[, 
                                                               1])
  el = tidyr::separate(data = el, col = 1, into = c("num", 
                                                    "denom"), sep = "___", remove = F)
  g = igraph::graph_from_edgelist(as.matrix(el[, 2:3]))
  igraph::E(g)$weight = el$Score
  s = igraph::strength(g)
  d = igraph::degree(g)
  nodeStrength = data.frame(Node = names(s), 
                            Str = as.numeric(s),
                            degreename = names(d),
                            degree = as.numeric(d),
                            meanDegree = scale(as.numeric(s)/as.numeric(d))) %>% 
    dplyr::arrange(dplyr::desc(meanDegree))
  nodeStrength
}
