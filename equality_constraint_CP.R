equality_constraint_CP = function(beta, GammaW, GammaM, delta_w, delta_m, 
                                  Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices) {
  
  NumGammaW = length(GammaW)
  NumGammaM = length(GammaM)
  NumGamma =  NumGammaW + NumGammaM
  
  # # omp_set_num_threads(4);
  
  probs = logitProb_latent_opp_set(beta, GammaW, GammaM, delta_w, delta_m, Xd, Zd, gw, gm, n, symmetric)
  CPW = probs$CPW
  CPM = probs$CPM
  
  sumVec = numeric(NumGamma)
  diffVec = numeric(NumGamma)
  
  # compute f(x_j,diamond) for all j women, where f(x_j,z_k) = w(x_j)*Px_jk * m(z_k)*Pz_kj
  for (j in 1:NumGammaW) {
    sum=0.0
    for (k in 1:NumGammaM) {
      sum = sum + CPW[j,k]*pmfM[k]*exp(gm) *  CPM[k,j]*pmfW[j]*exp(gw) * n 
      # sum = sum + CPW[j,k]*pmfM[k]*exp(gm) *  CPM[k,j]*pmfW[j]*exp(gw) * n / choices
    }
    sumVec[j]=sum
  }

  # compute f(diamond, z_j) for all j men, where f(z_j,x_k) = m(z_j)*Pz_jk * w(x_k)*Px_kj
  for (j in 1:NumGammaM) {
    sum=0.0
    for (k in 1:NumGammaW) {
      sum = sum + CPM[j,k]*pmfW[k]*exp(gw) * CPW[k,j]*pmfM[j]*exp(gm) * n
      # sum = sum + CPM[j,k]*pmfW[k]*exp(gw) * CPW[k,j]*pmfM[j]*exp(gm) * n /choices
    }
    sumVec[j+NumGammaW]=sum
  }

  # GammaW = f(x_j,diamond)/f(x_j,*), where f(x_j,*) = w(x_j)*Px_j0
  for (k in 1:NumGammaW) {
    # diffVec[k] = sumVec[k] - (pmfW[k] - (pmfW[k]*exp(gw)*CPW[k,1+NumGammaM])) 
    diffVec[k] = sumVec[k]/(pmfW[k]*exp(gw)*CPW[k,1+NumGammaM]) - (1-CPW[k,1+NumGammaM])/CPW[k,1+NumGammaM] # works
    # diffVec[k] = sumVec[k]/(pmfW[k]*exp(gw)*CPW[k,1+NumGammaM]) - GammaW[k] # works
  }
  # GammaM = f(diamond, z_j)/f(*,z_j), where f(*,z_j) = m(z_j)*Pz_j0
  for (k in 1:NumGammaM) {
    # diffVec[k+NumGammaW] = sumVec[k+NumGammaW] - (pmfM[k] - pmfM[k]*exp(gm)*CPM[k,1+NumGammaW])
    diffVec[k+NumGammaW] = sumVec[k+NumGammaW]/(pmfM[k]*exp(gm)*CPM[k,1+NumGammaW]) - (1-CPM[k,1+NumGammaW])/CPM[k,1+NumGammaW] # works
    # diffVec[k+NumGammaW] = sumVec[k+NumGammaW]/(pmfM[k]*exp(gm)*CPM[k,1+NumGammaW]) - GammaM[k] # works
  }
  
  # if (any(beta != 0))
  #    browser()
 
  
  # BLP contraction
  pmfj_est = matrix(0,nrow=1+NumGammaW, ncol=1+NumGammaM)
  # couples
  for (ii in 1:NumGammaW) {
    for (jj in 1:NumGammaM) {
      pmfj_est[ii, jj] = CPW[ii,jj]*pmfM[jj]*exp(gm) * CPM[jj,ii]*pmfW[ii]*exp(gw) * n / (exp(gw)+exp(gm))
      pmfj_est[ii, jj] = pmfj_est[ii, jj] * 2 # each couple consists of 2 people
    }
  }
  # single women
  for (ii in 1:NumGammaW) {
    pmfj_est[ii, 1+NumGammaM] = CPW[ii,1+NumGammaM]*pmfW[ii]*exp(gw) / (exp(gw)+exp(gm)) # / 2 / (exp(gw)+exp(gm))
  }
  # single men
  for (jj in 1:NumGammaM) {
    pmfj_est[1+NumGammaW, jj] = CPM[jj,1+NumGammaW]*pmfM[jj]*exp(gm) / (exp(gw)+exp(gm)) # / 2 / (exp(gw)+exp(gm))
  }
  
  # pmfj_est = pmfj_est * n # just get the raw number (sum(pmfj_est) = 2 because we set exp(gw) + exp(gm) = 2)
  
  matched_wm = log(sum(pmfj[1:NumGammaW, 1:NumGammaM])/sum(pmfj_est[1:NumGammaW, 1:NumGammaM]))
  # single_w = log(sum(pmfj[1:NumGammaW, 1+NumGammaM])/sum(pmfj_est[1:NumGammaW, 1+NumGammaM]))
  # single_m = log(sum(pmfj[1+NumGammaW, 1:NumGammaM])/sum(pmfj_est[1+NumGammaW, 1:NumGammaM]))
  # single_wm = log(sum(pmfj[1:NumGammaW, 1+NumGammaM])/sum(pmfj_est[1:NumGammaW, 1+NumGammaM]))
  
  # matched_wm = sum(pmfj[1:NumGammaW, 1:NumGammaM])-sum(pmfj_est[1:NumGammaW, 1:NumGammaM])
  # single_w = sum(pmfj[1:NumGammaW, 1+NumGammaM])-sum(pmfj_est[1:NumGammaW, 1+NumGammaM])
  # single_m = sum(pmfj[1+NumGammaW, 1:NumGammaM])-sum(pmfj_est[1+NumGammaW, 1:NumGammaM])
  # single_wm = sum(pmfj[1:NumGammaW, 1+NumGammaM], pmfj[1+NumGammaW, 1:NumGammaM])-
  #   sum(pmfj_est[1:NumGammaW, 1+NumGammaM], pmfj_est[1+NumGammaW, 1:NumGammaM])
    
  # if (delta_w != 0 || delta_m != 0)
  #    browser()
  
  # r_w = 1 # anything greater than 1.0e-6
  # while(abs(r_w) <= 1.0e-6) {
  #   r_w = log(matched_shares_w/predicted_matched_shares_w)
  #   matched_delta_w <<- matched_delta_w + r_w
  # }
  # r_m = 1 # anything greater than 1.0e-6
  # while(abs(r_m) <= 1.0e-6) {
  #   r_m = log(matched_shares_m/predicted_matched_shares_m)
  #   matched_delta_m <<- matched_delta_m + r_m
  # }
  
  # return(diffVec)
  # BLP contraction
  # return(c(diffVec, matched_wm, single_w, single_m))
  # return(c(diffVec, matched_wm, single_wm))
  return(c(diffVec, matched_wm))
  
}

