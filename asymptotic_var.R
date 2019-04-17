llik_paired = function(theta, NumBeta, NumGammaW, NumGammaM, j, k, pmfM, pmfW, Xd, Zd, gw, gm, n, symmetric) {
  
  beta = theta[1:NumBeta]
  GammaW = theta[NumBeta+(1:NumGammaW)]
  GammaM = theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
  delta_w = theta[NumBeta+NumGammaW+NumGammaM+1]
  delta_m = theta[NumBeta+NumGammaW+NumGammaM+2]
  
  probs = logitProb_latent_opp_set(beta, GammaW, GammaM, delta_w, delta_m, Xd, Zd, gw, gm, n, symmetric)
  CPW = probs$CPW
  CPM = probs$CPM
  
  llik = log(CPW[j,k]*pmfM[k]*exp(gm) * CPM[k,j]*pmfW[j]*exp(gw) * n / (exp(gw)+exp(gm)))
  return(llik)
}

llik_sw = function(theta, NumBeta, NumGammaW, NumGammaM, j, pmfW, Xd, Zd, gw, gm, n, symmetric) {
  
  beta = theta[1:NumBeta]
  GammaW = theta[NumBeta+(1:NumGammaW)]
  GammaM = theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
  delta_w = theta[NumBeta+NumGammaW+NumGammaM+1]
  delta_m = theta[NumBeta+NumGammaW+NumGammaM+2]
  
  probs = logitProb_latent_opp_set(beta, GammaW, GammaM, delta_w, delta_m, Xd, Zd, gw, gm, n, symmetric)
  CPW = probs$CPW
  CPM = probs$CPM
  
  llik = log(CPW[j,(1+NumGammaM)]*pmfW[j]*exp(gw) / (exp(gw)+exp(gm)))
  return(llik)
}

llik_sm = function(theta, NumBeta, NumGammaW, NumGammaM, k, pmfM, Xd, Zd, gw, gm, n, symmetric) {
  
  beta = theta[1:NumBeta]
  GammaW = theta[NumBeta+(1:NumGammaW)]
  GammaM = theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
  delta_w = theta[NumBeta+NumGammaW+NumGammaM+1]
  delta_m = theta[NumBeta+NumGammaW+NumGammaM+2]
  
  probs = logitProb_latent_opp_set(beta, GammaW, GammaM, delta_w, delta_m, Xd, Zd, gw, gm, n, symmetric)
  CPW = probs$CPW
  CPM = probs$CPM
  
  llik = log(CPM[k,(1+NumGammaW)]*pmfM[k]*exp(gm) / (exp(gw)+exp(gm))) 
  return(llik)
}

# not centered
outer_prod_gradient_CP = function(theta, NumBeta, NumGammaW, NumGammaM, Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric){
  
  outer_grad = matrix(0, nrow=length(theta), ncol=length(theta))
  
  # couples
  for (j in 1:NumGammaW) {
    for (k in 1:NumGammaM) {
      
      score = nl.grad(theta, llik_paired, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM, j=j, k=k, pmfM=pmfM, pmfW=pmfW, 
                      Xd=Xd, Zd=Zd, gw=gw, gm=gm, n=n, symmetric=symmetric)
      
      outer_grad = outer_grad + pmfj[j,k] * outer(score,score)
    }
  }
  
  # singles women (X each row)
  for (j in 1:NumGammaW) {
    
    score = nl.grad(theta, llik_sw, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM, j=j, pmfW=pmfW,  
                    Xd=Xd, Zd=Zd, gw=gw, gm=gm, n=n, symmetric=symmetric) 
    
    outer_grad = outer_grad + pmfj[j,(1+NumGammaM)] * outer(score,score)
  }
  
  # single men (Z each col)
  for (k in 1:NumGammaM) {
    
    score = nl.grad(theta, llik_sm, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM, k=k, pmfM=pmfM, 
                    Xd=Xd, Zd=Zd, gw=gw, gm=gm, n=n, symmetric=symmetric) 
    
    outer_grad = outer_grad + pmfj[(1+NumGammaW),k] * outer(score,score)
  }
  
  return (outer_grad)
  
}

# centered
outer_prod_gradient_CP_2 = function(theta, loglikfun, NumBeta, NumGammaW, NumGammaM, Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling){
  
  outer_grad = matrix(0, nrow=length(theta), ncol=length(theta))
  
  # couples
  for (j in 1:NumGammaW) {
    for (k in 1:NumGammaM) {
      score = nl.grad(theta, llik_paired, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM, j=j, k=k, pmfM=pmfM, pmfW=pmfW, 
                      Xd=Xd, Zd=Zd, gw=gw, gm=gm, n=n, symmetric=symmetric) - 
        nl.grad(theta, loglikfun, Xd=Xd, Zd=Zd, NumGammaW=NumGammaW, NumGammaM=NumGammaM,
                pmfW=pmfW, pmfM=pmfM,pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling, choices=1) # choice is not used
      
      outer_grad = outer_grad + pmfj[j,k] * outer(score,score)
    }
  }
  
  # singles women (X each row)
  for (j in 1:NumGammaW) {
    score = nl.grad(theta, llik_sw, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM, j=j, pmfW=pmfW,  
                    Xd=Xd, Zd=Zd, gw=gw, gm=gm, n=n, symmetric=symmetric) - 
      nl.grad(theta, loglikfun, Xd=Xd, Zd=Zd, NumGammaW=NumGammaW, NumGammaM=NumGammaM,
              pmfW=pmfW, pmfM=pmfM,pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling, choices=1) # choice is not used
    
    outer_grad = outer_grad + pmfj[j,(1+NumGammaM)] * outer(score,score)
  }
  
  # single men (Z each col)
  for (k in 1:NumGammaM) {
    score = nl.grad(theta, llik_sm, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM, k=k, pmfM=pmfM, 
                    Xd=Xd, Zd=Zd, gw=gw, gm=gm, n=n, symmetric=symmetric) - 
      nl.grad(theta, loglikfun, Xd=Xd, Zd=Zd, NumGammaW=NumGammaW, NumGammaM=NumGammaM,
              pmfW=pmfW, pmfM=pmfM,pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling, choices=1) # choice is not used
    
    outer_grad = outer_grad + pmfj[(1+NumGammaW),k] * outer(score,score)
  }
  
  return (outer_grad)
  
}



ave_hessian_CP = function(theta, NumBeta, NumGammaW, NumGammaM, Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric){
  
  ave_hess = matrix(0, nrow=length(theta), ncol=length(theta))
  
  # couples
  for (j in 1:NumGammaW) {
    for (k in 1:NumGammaM) {
      
      hess = hessian(llik_paired, theta, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM, j=j, k=k, pmfM=pmfM, pmfW=pmfW, 
                     Xd=Xd, Zd=Zd, gw=gw, gm=gm, n=n, symmetric=symmetric)
      
      ave_hess = ave_hess + pmfj[j,k] * hess
    }
  }
  
  # singles women (X each row)
  for (j in 1:NumGammaW) {
    
    hess = hessian(llik_sw, theta, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM, j=j, pmfW=pmfW,  
                   Xd=Xd, Zd=Zd, gw=gw, gm=gm, n=n, symmetric=symmetric)
    
    ave_hess = ave_hess + pmfj[j,(1+NumGammaM)] * hess
  }
  
  # single men (Z each col)
  for (k in 1:NumGammaM) {
    
    hess = hessian(llik_sm, theta, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM, k=k, pmfM=pmfM, 
                   Xd=Xd, Zd=Zd, gw=gw, gm=gm, n=n, symmetric=symmetric)
    
    ave_hess = ave_hess + pmfj[(1+NumGammaW),k] * hess
  }
  
  return (ave_hess)
  
}