simDataset = function( n = 100,
                       dims = 100,
                       baseNoise = 1e-11,
                       noiseCenter1 = -1, 
                       offset_c1 = -1, offset_c2 = -.5,
                       varDiag = .1,
                       num_shiftedComponents = 1,
                       num_noiseComponents = round(.75*(dims-num_shiftedComponents)),
                       size = 20,
                       avgReads1 = 1e7,
                       convertCounts = T){
 
  
  ##---------------------
  ##datset 1
  ##-------------------
  ## Class1
  m = rep(baseNoise,dims)
  beg = num_shiftedComponents+1
  end_ = (num_shiftedComponents+1)+(num_noiseComponents-1) 
  m[beg:end_] = noiseCenter1
  m[1:num_shiftedComponents] = offset_c1
  tt = clo(c(exp(m),1))
  message(log(tt[1]/tt[length(tt)])) ## expected
  s = diag(varDiag,dims)
  y1 = MASS::mvrnorm(n, mu = m, Sigma = s, tol = 1e-06, 
                     empirical = FALSE)
  y1 = alrInv(y1)
  colMeans(y1)
  y1[,1]/y1[,ncol(y1)]
  ### actual
  message(mean(y1[,1]/y1[,ncol(y1)]))
  hist((y1[,1]/y1[,ncol(y1)]))
  y1 = data.frame(y1)
  if(convertCounts){
    libSize = stats::rnbinom(n = nrow(y1), size = size, mu = avgReads1)
    y1 = round(sweep(y1,libSize,MARGIN = 1,FUN = "*"))
  }
  s1 = data.frame(Status = "S1",y1)
  
  ## Class 2
  m = rep(baseNoise,dims)
  beg = num_shiftedComponents+1
  end_ = (num_shiftedComponents+1)+(num_noiseComponents-1) 
  m[beg:end_] = noiseCenter1
  m[1:num_shiftedComponents] = offset_c2
  tt = clo(c(exp(m),1))
  message( tt[1]/tt[length(tt)]) ## expected
  s = diag(varDiag,dims)
  y1 = MASS::mvrnorm(n, mu = m, Sigma = s, tol = 1e-06, 
                     empirical = FALSE)
  
  y1 = alrInv(y1)
  colMeans(y1)
  y1[,1]/y1[,ncol(y1)]
  ### actual
  message(mean(y1[,1]/y1[,ncol(y1)]))
  hist((y1[,1]/y1[,ncol(y1)]))
  y1 = data.frame(y1)
  if(convertCounts){
    libSize = stats::rnbinom(n = nrow(y1), size = size, mu = avgReads1)
    y1 = round(sweep(y1,libSize,MARGIN = 1,FUN = "*"))
  }
  
  s2 = data.frame(Status = "S2",y1)
  
  ## MErge Data
  d1 = rbind(s1,s2)
}
