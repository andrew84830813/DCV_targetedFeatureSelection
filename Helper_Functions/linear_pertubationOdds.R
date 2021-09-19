#' Compute the Norm of a Closed Compositional Vector
#'
#'Compute the vector norm of a compositional vector in the simplex using Aitchinson's geometry
#'
#' @param namedVec a closed compositional vector
#'
#' @return the vector norm
#' @export
#'
vecNorm = function(namedVec){
  c = combinat::combn2(names(namedVec))
  el  = data.frame(c,ff = paste0(c[,1],"___",c[,2]))
  lrs = selEnergyPermR::getLogRatios(data.frame(t(namedVec)),ratioList = el)
  sqrt(sum(lrs^2)/length(lrs))
}



#' Distance Between Compositional Vectors
#'
#' Compute the distance between vectors using an all pairwise log ratio approach
#'
#' @param v1 closed compositional vector 1
#' @param v2 closed compositional vector 2
#'
#' @return the distance between vecotrs using all pariwise ratios
#' @export
#'
vecDist = function(v1,v2){
  c = combinat::combn2(names(v1))
  ff = paste0(c[,1],"___",c[,2])
  el  = data.frame(c,ff)
  
  v1 = selEnergyPermR::getLogRatios(data.frame((v1)),ratioList = el)
  
  v2 = selEnergyPermR::getLogRatios(data.frame((v2)),ratioList = el)
  
  sqrt(sum((v1-v2)^2)/length(v1))
}




#' Eigen decompositon of Composiitonal Dataset
#'
#' Performs eigen decompositon on compositional dataset after imputing zeroes and clr transform
#'
#' @param compMat a compositional data matrix/frame
#'
#' @return A list containing:\tabular{ll}{
#'    \code{clrmat} \tab clr transformed data  \cr
#'    \tab \cr
#'    \code{Coords} \tab Eigendecomposition coordinates \cr
#'     \tab \cr
#'    \code{eigenVector} \tab eigenvectors \cr
#'     \tab \cr
#'    \code{eigenValues} \tab eigenvalues \cr
#'     \tab \cr
#'    \code{percentExplained} \tab percent of variation explained by each component \cr
#' }
#' @export
#'
eigenDecomp.CLR <-
  function(compMat){
    X = as.matrix( compositions::clo(compMat) )    #<---- Data Matrix
    if(sum(X==0)>0){
      X = selEnergyPermR::fastImputeZeroes(X)
    }
    
    X.clr = data.frame(compositions::clr(X))
    X.clr.cov = stats::cov(X.clr)
    X.eigen = eigen(X.clr.cov)
    X.eigenvals = X.eigen$values
    X.evct = as.matrix(X.eigen$vectors)
    #compute score or pca coords
    coords = data.frame(as.matrix(X.clr)%*%X.evct)
    colnames(coords) = paste("Comp",1:ncol(coords),sep = "")
    perVar = (X.eigenvals/sum(X.eigenvals))*100
    
    X.evct = data.frame(X.evct,row.names = colnames(X))
    return(list(clrmat = X.clr,
                Coords = coords,
                eigenVector = X.evct,
                eigenValues = X.eigenvals,
                percentExplained = cumsum(perVar)))
  }






#' Linearly Shift a Point Cloud in the Simplex
#'
#'
#' Linearly shifts points cloud in the simplex along the directions of the centroid, most variation (eigenvector-1), or a user specified direction
#'
#' @importFrom foreach %dopar%
#'
#' @param tbl2 closed and imputed data matrix to be shifted within the simplex
#' @param pertubationDirection A user defined direction to shift point cloud; default = NULL
#' @param directionType The direction type of a standar shift in the simplex; either ("Centroid", "Eigen")
#' @param evComp if the directionType=="eigen" then which component to use for the linear shift; deafult = 1 for the first eigenvector
#' @param a A vector of \eqn{\alpha} power coeeficients to scale the linear shift; determines the distance of the shift
#' @param alpha_n Return the point cloud shift corresponding to a specific\eqn{\alpha} coefficient
#'
#' @return A list containing:\tabular{ll}{
#'    \code{combinedOutput} \tab xxx  \cr
#'    \tab \cr
#'    \code{linAdjustedData} \tab xxx \cr
#'     \tab \cr
#'    \code{allShiftData} \tab xxx \cr
#'     \tab \cr
#'    \code{shiftVec} \tab xxx \cr
#'     \tab \cr
#'    \code{meanDelta} \tab xxx \cr
#' }
#'
#' @export
#'
linearShift = function(tbl2,pertubationDirection = NULL,directionType = "centroid",evComp = 1,a = 1,alpha_n = 1){
  
  i = NULL
  ## centroids
  compMean.shifted = compositions::mean.acomp( compositions::acomp(tbl2) )
  
  ## Direction to mean
  if(is.null(pertubationDirection)){
    
    if(directionType == "centroid" || directionType== "Centroid" || directionType=="C"){
      pertb = compMean.shifted
    }else if(directionType =="EV" || directionType =="eigen"){
      c1 = eigenDecomp.CLR(tbl2)
      ev.c1  = c1$eigenVector
      c1.evDir = compositions::clrInv(ev.c1[,evComp])
      ## compute pertubation
      pertb = c1.evDir
    }
    
  }else{
    pertb = pertubationDirection
    names(pertb) = names(compMean.shifted)
  }
  
  ### -----------------------------------------*
  ## Alpha Range
  ### ------------------------------------------*
  shiftsPoints = foreach::foreach(i = 1:nrow(tbl2),.combine = rbind)%dopar%{
    c.df_pertb = data.frame()
    for(iA in a){
      ph = compositions::clo( t((as.matrix(tbl2[i,])))  * compositions::clo(pertb^iA) )
      c.df_pertb = rbind(c.df_pertb,ph)
    }
    data.frame(sampleNum = i,Alpha = a,c.df_pertb)
  }
  
  ## select point cloud relative to specific alpha
  pointAtAlpha_n = shiftsPoints[shiftsPoints$Alpha==alpha_n,-2:-1]
  
  
  
  simDat = rbind(
    data.frame(Type = "linearAdjusted",pointAtAlpha_n),
    data.frame(Type = "Base",tbl2))
  
  
  ## centroids
  compMean.new = compositions::mean.acomp( compositions::acomp(pointAtAlpha_n) )
  
  ## calculate true distance between centroids
  
  ds = vecDist(data.frame(t(compMean.new)),data.frame(t(compMean.shifted))) #sqrt(sum( ( as.numeric(clr(compMean.new)) - as.numeric(clr(compMean.shifted)) )^2   ))
  
  return(list(combinedOutput = simDat,
              linAdjustedData = pointAtAlpha_n,
              allShiftData = shiftsPoints,
              shiftVec = ( as.numeric( compositions::clr(compMean.new)) - as.numeric( compositions::clr(compMean.shifted)) )^2 ,
              meanDelta = c(alpha = a,distance=ds))
  )
}




compute_oddsAUC = function(odds_table,posClass,negClass){
  ## compute areas
  or.df1 = odds_table %>% 
    group_by(alpha,Part,Status) %>% 
    summarise_all(.funs = mean)
  
  part = unique(odds_table$Part)
  
  ## 
  increase.abund = or.df1 %>% 
    filter(alpha>=0)
  auc.increase_posClass = MESS::auc(x = increase.abund$alpha[increase.abund$Status==posClass],
                                    y = increase.abund$logodds_pos[increase.abund$Status==posClass]
  )
  auc.increase_negClass = MESS::auc(x = increase.abund$alpha[increase.abund$Status==negClass],
                                    y = increase.abund$logodds_pos[increase.abund$Status==negClass]
  )
  
  decrease.abund = or.df1 %>% 
    filter(alpha<=0)
  auc.decrease_posClass = MESS::auc(x = decrease.abund$alpha[decrease.abund$Status==posClass],
                                    y = decrease.abund$logodds_pos[decrease.abund$Status==posClass]
  )
  auc.decrease_negClass = MESS::auc(x = decrease.abund$alpha[decrease.abund$Status==negClass],
                                    y = decrease.abund$logodds_pos[decrease.abund$Status==negClass]
  )
  
  data.frame(part,auc.increase_posClass,auc.increase_negClass,auc.decrease_posClass,auc.decrease_negClass)
}

