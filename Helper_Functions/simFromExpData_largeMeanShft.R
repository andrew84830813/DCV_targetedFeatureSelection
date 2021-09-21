


simFromExpData.largeMeanShft <-
  function( raMatrix ,n1 = 30 , n2 = 30,perFixedFeatures = 0.95,featureShiftPercent = 1.25,seed = 08272008){
    mu = NULL
    var = NULL
    set.seed(seed)
    y.alr =  compositions::alr( raMatrix )
    alr.cov = stats::cov(y.alr)
    var.alr = data.frame(ratio = names(diag(alr.cov)),var = diag(alr.cov)) %>%
      dplyr::arrange(-var)
    y.alr = subset(data.frame(y.alr), select = var.alr$ratio)
    alr.mean = colMeans(y.alr)
    alr.cov = stats::cov(y.alr)
    cv = mean(abs( stats::sd(var.alr$var)/alr.mean))
    
    
    alr.sim = MASS::mvrnorm(n = n1,mu = alr.mean,Sigma = alr.cov)
    alr.mean2 = colMeans(alr.sim)
    alr.cov = stats::cov(alr.sim)
    alr.ra =  compositions::alrInv(alr.sim); colnames(alr.ra) = colnames(raMatrix)
    s1 = data.frame(Status = "S1",alr.ra)
    
    ####
    nfeats = length(alr.mean)
    hold = round(length(alr.mean2)*perFixedFeatures)
    alr.mean2[(hold+1):nfeats] = featureShiftPercent * alr.mean[(hold+1):nfeats]
    alr.sim = MASS::mvrnorm(n = n2,mu = alr.mean2,Sigma = alr.cov)
    alr.ra =  compositions::alrInv(alr.sim); colnames(alr.ra) = colnames(raMatrix)
    s2 = data.frame(Status = "S2",alr.ra)
    ###
    
    
    rbind(s1,s2)
  }
