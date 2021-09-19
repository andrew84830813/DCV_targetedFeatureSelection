s1 = function (n1 = 30, n2 = 30, dms_ = 75, seed,alphaScale = 2) 
{
  a = rep(1, dms_)
  a = a
  expectedScale = dms_
  scale = sum(a)
  scaleRatio = scale/(dms_)
  set.seed(seed)
  s1 = selEnergyPermR::sampleDirichlet(n1 = n1, dims = dms_, sampleName = "S1", 
                       a1 = a * 3)$Sample
  s2 = selEnergyPermR::sampleDirichlet(n1 = n2, dims = dms_, sampleName = "S2", 
                       a1 = a * alphaScale)$Sample
  df = rbind(s1, s2)
  return(df)
}




s3 = function (n1 = 30, n2 = 30, dms_ = 75, seed,muFactor = 1,perFeatures = .25) 
{
  set.seed(seed)
  mu_ = rep(0, dms_)
  t25 = round(dms_ * perFeatures)
  mu2 = rep(0, dms_)
  mu2[1:t25] = muFactor
  sigma_ = diag(dms_)
  # diag(sigma_[-1, ]) <- rep(0.2, ncol(sigma_) - 1)
  # sigma_ = t(sigma_)
  # diag(sigma_[-1, ]) <- rep(0.2, ncol(sigma_) - 1)
  sigma_ = Matrix::nearPD(sigma_)$mat
  U = matrix(stats::runif(dms_ * dms_, 0, 32/(dms_^2)), nrow = dms_)
  U_ = sigma_ + U
  U_ = Matrix::nearPD(U_)$mat
  eg = min(eigen(as.matrix(sigma_))$values)
  sig = min(eigen(as.matrix(U_))$values)
  dd = min(eg, sig) + 0.05
  sig1 = Matrix::nearPD(sigma_ + dd * diag(dms_))$mat
  sig2 = Matrix::nearPD(sigma_ + U + dd * diag(dms_))$mat
  s1 = sampleAddLogisticNormal(n = n1, dims_ = dms_, mu = mu_, 
                               sigma = sig1, sigmaScale = 1, sampleName = "S1")
  s2 = sampleAddLogisticNormal(n = n2, dims_ = dms_, mu = mu2, 
                               sigma = sig1, sigmaScale = 1, sampleName = "S2")
  df = rbind(s1, s2)
}







