#' Reconstruct the matching based on the estimated parameters
#' 
#' Calculate the joint density of the observed features of the pairs between agents
#' using the estimated parameters, given the feature matrices of the two sides. The 
#' observed joint density is also computed.  
#' 
#' @param ff formula; an \code{\link{formula}} object, of the form \code{
#' ~ <model terms>}. For the details on the possible \code{<model terms>}, see
#' \code{\link{rpm-terms}}.
#' @param theta The estimated parameters.
#' @param mu The observed matching matrix, where 1 represents a pairing, 0 otherwise. 
#' Each row is a woman, each column is a man. The order of the rows (columns)
#' needs to be the same as in \code{X} (\code{Z}). 
#' @param X Feature matrix for women. Each row is a woman, each column is a feature. 
#' The number of column is assumed to be the same as in \code{Zdata}. 
#' @param Z Feature matrix for men. Each row is a man, each column is a feature. 
#' The number of column is assumed to be the same as in \code{Xdata}.
#' @param symmetric If set to \code{TRUE}, the same utility coefficients are to be used 
#' for both the women's and men's sides; otherwise, separate estimates need to be provided for each side. 
#' @return This function returns a list consisting of the following elements: 
#' \item{pmfj_est}{The reconstructed joint density matrix based on the estimated parameters.}
#' \item{pmfj_obs}{The observed joint density matrix based on the data \code{X} and \code{Z}.}
#' @seealso fitrpm_R_CP
#' @references Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @keywords models
#' @examples
#' 
check_CP_latent = function(ff, theta, mu, X, Z, X_w, Z_w, pair_w, symmetric, choices){

  # num_single_women = nrow(modelmat$X) - sum(rowSums(mu)) 
  # num_single_men = nrow(modelmat$Z) - sum(colSums(mu)) 
  # n=nrow(modelmat$X)+nrow(modelmat$Z)
  n = nrow(X) + nrow(Z)

  Xdata <- cbind(1, X)
  Zdata <- cbind(1, Z)
  colnames(Xdata)[1] <- "Int"
  colnames(Zdata)[1] <- "Int"
  
  model.terms <- rownames(attr(terms.formula(ff), "factors"))
  temp <- strsplit(model.terms, "[(]")
  model.terms.names <- unlist(lapply(temp, `[[`, 1))
  temp2 <- gsub(")", "", gsub("\"", "", unlist(temp)))
  temp2 <- gsub("[ \t\n\r\f\v]", "", temp2)
  model.terms.coef.names <- temp2[-which(temp2 %in% model.terms.names)]
  model.terms.coef.names <- strsplit(model.terms.coef.names, ",")

  # 2) get the subset of the variables that are relevant according to the formula
  # Define the variables used in the model (and hence the unique classes of partners)
  # This is typically a subset of the available variables
  model_vars <- c("Int", unlist(unique(lapply(model.terms.coef.names, '[[', 1))))

  Xu <- unique(Xdata[,model_vars])
  Xu <- Xu[do.call(order, as.data.frame(Xu)),]
  Zu <- unique(Zdata[,model_vars])
  Zu <- Zu[do.call(order, as.data.frame(Zu)),]
  
  Xtype <- rep(NA,nrow(Xdata))
  for(i in 1:nrow(Xu)){
    Xtype[apply(Xdata[,model_vars], 1, function(x) identical(x, Xu[i,]))] <- i
  }
  # Ztype: group membership for men (one for each man in the pop)
  Ztype <- rep(NA,nrow(Zdata))
  for(i in 1:nrow(Zu)){
    Ztype[apply(Zdata[,model_vars], 1, function(x) identical(x, Zu[i,]))] <- i
  }
  
  pmfW = wtd.table(Xtype, weights=X_w)
  pmfW = pmfW/sum(pmfW)
  pmfM = wtd.table(Ztype, weights=Z_w)
  pmfM = pmfM/sum(pmfM)
  
  num_Xu = nrow(Xu)
  num_Zu = nrow(Zu)
  
  
  ########## modified to work with many-t0-many ###########
  ####### only works for "INDIV" for now #############
  # this may be depend on the application
  # the code below assumes mu is a diagonal matrix and there are no singles
  ####### only works for "INDIV" for now #############
  match_w = apply(mu, 1, function(x) which(x>=1))
  match_w = lapply(match_w, function(x){Ztype[x, drop=FALSE]})
  
  # match_w = lapply(match_w, function(x){c(x,rep((num_Zu+1), choices-length(x)))})
  # w_table = table(rep(Xtype, each=choices), unlist(match_w))
  match_w = unlist(match_w)
  df = data.frame(type = Xtype, m_type=match_w, freq=mu@x, weight=pair_w)
  df.exp = df[rep(row.names(df), df$freq), 1:4]
  dfs = data.frame(type = Xtype, m_type=rep(num_Zu+1, length(Xtype)), freq=choices - mu@x, weight=pair_w)
  dfs = dfs[dfs$freq>0,]
  dfs.exp = dfs[rep(row.names(dfs), dfs$freq), 1:4]
  # w_table = table(c(df.exp$type, dfs.exp$type), c(df.exp$m_type, dfs.exp$m_type))
  w_table = wtd.table(c(df.exp$type, dfs.exp$type), c(df.exp$m_type, dfs.exp$m_type), weights=c(df.exp$weight, dfs.exp$weight))
  
  w_table = rbind(w_table, 0)
  
  match_m = apply(mu, 2, function(x) which(x>=1))
  match_m = lapply(match_m, function(x){Xtype[x, drop=FALSE]})
  
  # match_m = lapply(match_m, function(x){c(x,rep((num_Xu+1), choices-length(x)))})
  # m_table = table(rep(Ztype, each=choices), unlist(match_m))
  match_m = unlist(match_m)
  df = data.frame(type = Ztype, m_type=match_m, freq=mu@x, weight=pair_w)
  df.exp = df[rep(row.names(df), df$freq), 1:4]
  dfs = data.frame(type = Ztype, m_type=rep(num_Xu+1, length(Ztype)), freq=choices - mu@x, weight=pair_w)
  dfs = dfs[dfs$freq>0,]
  dfs.exp = dfs[rep(row.names(dfs), dfs$freq), 1:4]
  # m_table = table(c(df.exp$type, dfs.exp$type), c(df.exp$m_type, dfs.exp$m_type))
  m_table = wtd.table(c(df.exp$type, dfs.exp$type), c(df.exp$m_type, dfs.exp$m_type), weights=c(df.exp$weight, dfs.exp$weight))
  
  m_table = rbind(m_table, 0)
  
  pmfj = w_table + t(m_table) # just raw number (pairs are doubled by adding up w and m)
  pmfj = pmfj/sum(pmfj)
  
  
  modelmat <- rpm.model.matrix(model.terms.names, model.terms.coef.names, Xu, Zu)

  # assume same # of explanatory variables for both men's and women's side
  NumGammaW <- num_Xu
  NumGammaM <- num_Zu
  NumGamma <- NumGammaW + NumGammaM
  if (symmetric) NumBeta <- dim(modelmat$X)[3] else NumBeta <- dim(modelmat$X)[3]+dim(modelmat$Z)[3]
  beta <- theta[1:NumBeta]
  GammaW <- theta[NumBeta+(1:NumGammaW)]
  GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
  delta_w = theta[NumBeta+NumGammaW+NumGammaM+1]
  delta_m = theta[NumBeta+NumGammaW+NumGammaM+2]
  
  # get the proportion of men and women
  gw = log(nrow(Xdata)/n*2) 
  gm = log(nrow(Zdata)/n *2) 
  # gw = gm = 0
  
  CP_result=logitProb_latent_opp_set(beta, GammaW, GammaM, delta_w, delta_m, modelmat$X, modelmat$Z, gw, gm, n, symmetric)
  
  CP_result$CPW
  CP_result$CPM
  
  # check joint probabilities
  # couples
  pmfj_est = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu)
  for (ii in 1:num_Xu) {
    for (jj in 1:num_Zu) {
      pmfj_est[ii, jj] = CP_result$CPW[ii,jj]*pmfM[jj]*exp(gm) * CP_result$CPM[jj,ii]*pmfW[ii]*exp(gw) * n / (exp(gw)+exp(gm))
      pmfj_est[ii, jj] = pmfj_est[ii, jj] * 2 # each couple consists of 2 people
    }
  }
  # single women
  for (ii in 1:num_Xu) {
    pmfj_est[ii, 1+num_Zu] = CP_result$CPW[ii,1+num_Zu]*pmfW[ii]*exp(gw) / (exp(gw)+exp(gm))
  }
  # single men
  for (jj in 1:num_Zu) {
    pmfj_est[1+num_Xu, jj] = CP_result$CPM[jj,1+num_Xu]*pmfM[jj]*exp(gm) / (exp(gw)+exp(gm))
  }
  
  # pmfj_est = pmfj_est * n # just get the raw number (sum(pmfj_est) = 2 because we set exp(gw) + exp(gm) = 2)
  
  
  # # debug
  # # sum of type 1 females ( ~= exp(gw) * pmfW[1] )
  # print(paste0("sum of type 1 females (should be ", exp(gw)*pmfW[1], " :"))
  # print(sum(pmfj_est[1,]))
  # 
  # # sum of type 2 females ( ~= exp(gw)  * pmfW[2] )
  # print(paste0("sum of type 2 females (should be ", exp(gw)*pmfW[2], " :"))
  # print(sum(pmfj_est[2,]))
  # 
  # # sum of all females ( ~= exp(gw) )
  # print(paste0("sum of all females (should be ", exp(gw), " :"))
  # print(sum(pmfj_est[1:2,]))
  # 
  # # sum of type 1 males ( ~= exp(gm) * pmfM[1] )
  # print(paste0("sum of type 1 males (should be ", exp(gm)*pmfM[1], " :"))
  # print(sum(pmfj_est[,1]))
  # 
  # # sum of type 2 males ( ~= exp(gm) * pmfM[2] )
  # print(paste0("sum of type 2 males (should be ", exp(gm)*pmfM[2], " :"))
  # print(sum(pmfj_est[,2]))
  # 
  # # sum of all males ( ~= exp(gm) )
  # print(paste0("sum of all males (should be ", exp(gm), " :"))
  # print(sum(pmfj_est[,1:2]))
  
  # # debug
  # # sum of population ( ~= 1 )
  # print("sum of population (should be 1")
  # print(sum(pmfj_est))

  # # compare estimated joint probabilities with truth
  # print("estimated joint probabilities")
  # print(pmfj_est)
  # print("observed joint probabilities")
  # print(pmfj)
  
  # # print likelihood value
  # loglik = loglikelihood_CP(beta, GammaW, GammaM, modelmat$X, modelmat$Z, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling="INDIV")
  # logliknull = loglikelihood_null_CP(GammaW, GammaM, modelmat$X, modelmat$Z, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling="INDIV")
  # logliknull_gamma = loglikelihood_null_CP_gamma(GammaW, GammaM, modelmat$X, modelmat$Z, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling="INDIV")
  # print(paste0("loglik = ", loglik))
  # print(paste0("logliknull = ", logliknull))
  # # print(paste0("logliknull_gamma = ", logliknull_gamma))
  
  # # chi-squared statistics
  # if (symmetric) {
  #   df = NumBeta
  # } else {
  #   df = NumBeta # same number of explanatory variables for both sides
  # }
  # # p.val1 <- 1 - pchisq(-2*(loglik-logliknull), df=df) 
  # # print(p.val1)
  # p.val2 <- pchisq(-2*(loglik-logliknull), df=df, lower.tail=FALSE) 
  # print(paste0("chi-squared p-val = ", p.val2))
  
  # # debug
  # diffVecnull_gamma = equality_constraint_null_CP_gamma(GammaW, GammaM, modelmat$X, modelmat$Z, pmfW, pmfM, gw, gm, n, symmetric, sampling="INDIV")
  # diffVecnull = equality_constraint_null_CP(GammaW, GammaM, modelmat$X, modelmat$Z, pmfW, pmfM, gw, gm, n, symmetric, sampling="INDIV")
  # print("diffVecnull_gamma ")
  # print(diffVecnull_gamma)
  # print("diffVecnull")
  # print(diffVecnull)
  
  return(list(pmfj_est=pmfj_est, pmfj_obs=pmfj))
  
}

find_best_theta = function(beta_med, beta_sim, numBeta, numGamma) {
  best_idx = 0
  best_diff = .Machine$double.xmax
  for (ii in 1:nrow(beta_sim)) {
    diffval = sum((beta_sim[ii,1:numBeta] - beta_med)^2)
    if (diffval < best_diff) {
      best_diff = diffval
      best_idx = ii
    }
  }
  
  print(paste0("best theta [", best_idx, "] = "))
  print(beta_sim[best_idx,1:(numBeta+numGamma)])
  
  return(beta_sim[best_idx,1:(numBeta+numGamma)])
}