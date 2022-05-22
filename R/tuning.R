cv_score_multigraph <- function(cv_data, PreS){
  d = dim(cv_data)[1]
  n = dim(cv_data)[2]
  p = dim(cv_data)[3]
  q = dim(cv_data)[4]
  data12 = array(aperm(cv_data, c(2,3,4,1)), dim = c(n*p, q, d))
  data12 = aperm(data12, c(3,1,2))
  score = 0
  for (t in 1:d){
    tempdata = data12[t,,]
    Sample_covS = cov(tempdata)
    score = score + compute_gaussian_loglike(Sample_covS, PreS[t,,])
  }
  return(score)
}

compute_gaussian_loglike <- function(S_sample, PreSt){
  ret = eigen(PreSt)
  V = ret$vectors
  lm = ret$values
  if (min(lm)<0){
    cat('min less than 0', min(lm), 'max is', max(lm), '\n')
    lm = lm - 2*min(lm)
    PreSt = V%*%diag(lm)%*%t(V)
  }
  loglike = log((det(PreSt))) - sum(diag(S_sample%*%PreSt))
  return(loglike)
}

auto_lambda <- function(yall, d, n, p, q, mode = "mean"){
  a1 = 0
  a2 = 0
  for (t in 1:d){
    yday = yall[t,,,]
    data12 = matrix(yday, n*p, q)
    data13 = aperm(yday, c(1,3,2))
    data13 = matrix(data13, n*q, p)
    covT = cov(data13)
    covT = covT/sum(diag(covT))*p
    covS = cov(data12)
    if (mode=='mean'){
      a1 = max(a1,mean(eigen(covT)$values))
      a2 = max(a2, mean(abs(diag(covS))))      
    }
    else{
      a1 = max(a1,max(eigen(covT)$values))
      a2 = max(a2, max(abs(diag(covS))))       
    }
  }
  ld = a1*a2*sqrt((d+log(q))/n/p)
  return(ld)
}

sample_est <- function(yall, d, n, p, q){
  a1 = 0
  a2 = 0
  samplerho = array(0, dim=c(d, q, q))
  samplePreS = array(0, dim=c(d, q, q))
  sampleCovT = array(0, dim=c(d, p, p))
  for (t in 1:d){
    yday = yall[t,,,]
    data12 = matrix(yday, n*p, q)
    data13 = aperm(yday, c(1,3,2))
    data13 = matrix(data13, n*q, p)
    covT = cov(data13)
    covT = covT/sum(diag(covT))*p
    covS = cov(data12)
    preS = solve(covS)
    sampleCovT[t,,] = covT
    samplePreS[t,,] = preS
    for (i in 1:q){
      for (j in 1:q){
        samplerho[t, i, j] = samplePreS[t, i, j]/sqrt(samplePreS[t, i, i]*samplePreS[t, j, j])
      }
    }
  }
  return(list(rho=samplerho, ct = sampleCovT, ps = samplePreS))
}