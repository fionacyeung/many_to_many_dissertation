#' Estimate the parameters of a Revealed Preference Matchings Model
#' 
#' \code{\link{fitrpm}} estimates parameters for men and women of certain
#' characteristics (or shared characteristics) of people of the opposite sex.
#' It does this using an approximate likelihood based on ideas from Menzel (2015).
#' 
#' @aliases rpm.object
#' @param ff formula; an \code{\link{formula}} object, of the form \code{
#' ~ <model terms>}. For the details on the possible \code{<model terms>}, see
#' \code{\link{rpm-terms}}.
#' @param mu The observed matching matrix, where 1 represents a pairing, 0 otherwise. 
#' Each row is a woman, each column is a man. The order of the rows (columns)
#' needs to be the same as in \code{Xdata} (\code{Zdata}). 
#' @param Xdata Feature matrix for women. Each row is a woman, each column is a feature. 
#' The number of column is assumed to be the same as in \code{Zdata}. 
#' @param Zdata Feature matrix for men. Each row is a man, each column is a feature. 
#' The number of column is assumed to be the same as in \code{Xdata}.
#' @param theta_0 vector; initial parameter values to initiate the search for the MLE. 
#' The length of the vector corresponds to the number of explanatory variables (number 
#' of columns in \code{Xdata}). Default vaule is a vector where the utility coefficients
#' are set to 0 and the expected utility of the opportunity set (for each distinct type of 
#' women and men) set to 1.
#' @param control A list of control parameters for algorithm tuning. Constructed using
#' \code{\link{control.fitmle}}. The \code{symmetric} parameter, if set to \code{TRUE}, indicates 
#' that the same utility coefficients are to be used for both the women's and men's sides; 
#' otherwise, they will be estimated separately for each side. The \code{sampling_protocol} parameter 
#' take the following values: \code{"INDV"} (default) (individuals are sampled, data contains both
#' singles and couples), \code{"COUPLE"} (only couples are included in the data), 
#' \code{"HOUSEHOLD"} (households are sampled, each household can be a single or a couple)
#' @return \code{\link{fitrpm}} returns an object of class \code{\link{rpm.object}}
#' that is a list consisting of the following elements: 
#' \item{coef}{The maximum likelihood estimate of \eqn{\theta}, the vector of
#' coefficients for the model parameters. This includes the model \eqn{\beta} and the model \eqn{\Gamma}.}
#' \item{loglik}{The value of the maximized log-likelihood.}
#' \item{exitflag}{integer value with the status of the optimization (4 is success).}
#' \item{call}{the call that was made to nloptr.}
#' \item{x0}{vector with starting values for the optimization.}
#' \item{message}{more informative message with the status of the optimization.}
#' \item{iterations}{number of iterations that were executed.}
#' \item{objective}{value if the objective function in the solution.}
#' \item{solution}{optimal value of the controls.}
#' \item{version}{version of NLopt that was used.}
#' \item{covar}{Approximate covriance matrix of the estimates.}
#' \item{eq}{Values from the equality constraints. Larger values indicate non-convergence.}
#' @seealso control.fitrpm, summary.rpm, print.rpm
#' @references Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @keywords models
#' @examples

library(abind)
library(nloptr)
library(numDeriv)
library(MASS)

fitrpm_R_CP <- function(formula, mu, Xdata, Zdata, theta_0=NULL, choices, control){

    symmetric = control[["symmetric"]]
    sampling = control[["sampling_protocol"]]
    
    # number of pairs in the data
    K=nrow(Xdata)
    
    # get the proportion of men and women
    n = nrow(Xdata) + nrow(Zdata)
    gw = log(nrow(Xdata)/n*2) # to ensure exp(gw)+exp(gm) = 2
    gm = log(nrow(Zdata)/n*2)
    
    ##############???????????
    # gw = gm = 0
    ######################
   
    # 1) parse the formula
    # intercept is always added as the first column
    Xdata <- cbind(1, Xdata)
    Zdata <- cbind(1, Zdata)
    colnames(Xdata)[1] <- "Int"
    colnames(Zdata)[1] <- "Int"
    
    model.terms <- rownames(attr(terms.formula(formula), "factors"))
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
    
    
    # 3) Compute the marginal distributions of women's and men's types
    Xu <- unique(Xdata[,model_vars])
    Xu <- Xu[do.call(order, as.data.frame(Xu)),]
    Zu <- unique(Zdata[,model_vars])
    Zu <- Zu[do.call(order, as.data.frame(Zu)),]
    
    
    # 4) Create joint PMF
    # Xtype: group membership for women (one for each woman in the pop)
    Xtype <- rep(NA,nrow(Xdata))
    for(i in 1:nrow(Xu)){
        Xtype[apply(Xdata[,model_vars], 1, function(x) identical(x, Xu[i,]))] <- i
    }
    # Ztype: group membership for men (one for each man in the pop)
    Ztype <- rep(NA,nrow(Zdata))
    for(i in 1:nrow(Zu)){
        Ztype[apply(Zdata[,model_vars], 1, function(x) identical(x, Zu[i,]))] <- i
    }
    
    # order the data by pair
    Xtype_paired = Xtype[unlist(apply(mu, 2, function(x) which(x>0)))]
    ##### modified to work for many-to-many
    # Ztype_paired = Ztype[as.logical(colSums(mu))]
    Ztype_paired = Ztype[unlist(apply(mu, 1, function(x) which(x>0)))] 
    
    # Xtype_single = table(Xtype[!rowSums(mu)])
    # Ztype_single = table(Ztype[!colSums(mu)])
    Xtype_single = table(factor(Xtype[!rowSums(mu)], 1:nrow(Xu))) # account for missing types
    Ztype_single = table(factor(Ztype[!colSums(mu)], 1:nrow(Zu))) # account for missing types
    
    pmfW=table(Xtype)/nrow(mu)
    pmfM=table(Ztype)/ncol(mu)
    
    num_Xu = nrow(Xu)
    num_Zu = nrow(Zu)

    ########## modified to work with many-t0-many ###########
    ####### only works for "INDIV" for now #############
    match_w = apply(mu, 1, function(x) which(x>=1))
    match_w = lapply(match_w, function(x){Ztype[x, drop=FALSE]})
    
    # match_w = lapply(match_w, function(x){c(x,rep((num_Zu+1), choices-length(x)))})
    # w_table = table(rep(Xtype, each=choices), unlist(match_w))
    match_w = unlist(match_w)
    df = data.frame(type = Xtype, m_type=match_w, freq=mu@x)
    df.exp = df[rep(row.names(df), df$freq), 1:3]
    dfs = data.frame(type = Xtype, m_type=rep(num_Zu+1, length(Xtype)), freq=choices - mu@x)
    dfs = dfs[dfs$freq>0,]
    dfs.exp = dfs[rep(row.names(dfs), dfs$freq), 1:3]
    w_table = table(c(df.exp$type, dfs.exp$type), c(df.exp$m_type, dfs.exp$m_type))
    
    w_table = rbind(w_table, 0)
    
    match_m = apply(mu, 2, function(x) which(x>=1))
    match_m = lapply(match_m, function(x){Xtype[x, drop=FALSE]})
    
    # match_m = lapply(match_m, function(x){c(x,rep((num_Xu+1), choices-length(x)))})
    # m_table = table(rep(Ztype, each=choices), unlist(match_m))
    match_m = unlist(match_m)
    df = data.frame(type = Ztype, m_type=match_m, freq=mu@x)
    df.exp = df[rep(row.names(df), df$freq), 1:3]
    dfs = data.frame(type = Ztype, m_type=rep(num_Xu+1, length(Ztype)), freq=choices - mu@x)
    dfs = dfs[dfs$freq>0,]
    dfs.exp = dfs[rep(row.names(dfs), dfs$freq), 1:3]
    m_table = table(c(df.exp$type, dfs.exp$type), c(df.exp$m_type, dfs.exp$m_type))
    
    m_table = rbind(m_table, 0)
    
    pmfj = w_table + t(m_table) # just raw number (pairs are doubled by adding up w and m)
    pmfj = pmfj/sum(pmfj)
    # pmfj = pmfj/n
    
    
    # make multinomial coeff from choice sequence
    c_seq_w = t(sapply(match_w, function(x){table(factor(x, 1:(num_Zu+1)))}))
    c_seq_m = t(sapply(match_m, function(x){table(factor(x, 1:(num_Xu+1)))}))
    mult_coeff_w = factorial(choices)/apply(factorial(c_seq_w),1,prod)
    mult_coeff_m = factorial(choices)/apply(factorial(c_seq_m),1,prod)
    
    ##################################################
    
    # if (sampling == "COUPLE") { 
    #   
    #   pmfj = matrix(0,nrow=num_Xu, ncol=num_Zu) # women (X) indexed by row, men (Z) indexed by column
    #   pmfj = unclass(table(Xtype_paired,Ztype_paired))
    #   pmfj = pmfj / nrow(Xdata)
    #   
    #   
    # } else if (sampling == "HOUSEHOLD") {
    #   
    #   pmfj = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
    #   pmfj[1:num_Xu,1:num_Zu] = unclass(table(Xtype_paired,Ztype_paired))
    #   
    #   if (length(Xtype_single) > 0) {
    #     pmfj[1:num_Xu,1+num_Zu] = Xtype_single
    #   } 
    #   if (length(Ztype_single) > 0) {
    #     pmfj[1+num_Xu,1:num_Zu] = Ztype_single
    #   }
    #   
    #   pmfj = pmfj / nrow(Xdata)
    #   
    #   
    # } else { # assume "INDIV"
    #   
    #   pmfj = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
    #   pmfj[1:num_Xu,1:num_Zu] = unclass(table(Xtype_paired,Ztype_paired)) *2
    #   
    #   if (length(Xtype_single) > 0) {
    #     pmfj[1:num_Xu,1+num_Zu] = Xtype_single
    #   } 
    #   if (length(Ztype_single) > 0) {
    #     pmfj[1+num_Xu,1:num_Zu] = Ztype_single
    #   }
    #   browser()
    #   ##### modified to work for many-to-many
    #   # pmfj = pmfj / (nrow(Xdata) + nrow(Zdata)) 
    #   pmfj = pmfj/(sum(pmfj))
    #   
    # }
   #############################################

    
    # 5) create model matrix
    modelmat <- rpm.model.matrix(model.terms.names, model.terms.coef.names, Xu, Zu)
    
    X <- modelmat$X
    Z <- modelmat$Z
    
  

    NumBetaW <- dim(X)[3]
    NumBetaM <- if(symmetric)0 else dim(Z)[3]
    NumBeta <- NumBetaW + NumBetaM
    NumGammaW <- num_Xu
    NumGammaM <- num_Zu
    NumGamma <- NumGammaW + NumGammaM
    
    
    # 6) Set theta_0 to a good starting value to save time
    if(is.null(theta_0)){
        theta_0 <- rep(c(0,1,0),c(NumBeta,NumGamma,2)) # 2 for delta
    }
    if(symmetric) {
      nstr = gsub("b1","", modelmat$Xnames)
    }else {
      nstr = c(modelmat$Xnames, modelmat$Znames)
    }
    names(theta_0) <- c(nstr, paste("ExpUtil.b1.",(1:NumGammaW),sep=""), paste("ExpUtil.b2.",(1:NumGammaM),sep=""), 
                        "delta_w", "delta_m")
    
    # core algorithm begins
    
    betaW <- theta_0[1:NumBetaW]
    betaM <- if(symmetric)NULL else theta_0[NumBetaW + (1:NumBetaM)]
    
    loglikfun <- function(theta, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices){
      if (symmetric) NumBeta <- dim(Xd)[3] else NumBeta <- dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      delta_w = theta[NumBeta+NumGammaW+NumGammaM+1]
      delta_m = theta[NumBeta+NumGammaW+NumGammaM+2]
      -loglikelihood_CP(beta, GammaW, GammaM, delta_w, delta_m, Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices)
    } 
    eqfun <- function(theta, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices){
      if (symmetric) NumBeta <- dim(Xd)[3] else NumBeta <- dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      delta_w = theta[NumBeta+NumGammaW+NumGammaM+1]
      delta_m = theta[NumBeta+NumGammaW+NumGammaM+2]
      equality_constraint_CP(beta, GammaW, GammaM, delta_w, delta_m, Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices)
    }
    gloglikfun <- function(theta, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices){
      nl.grad(theta, loglikfun,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
              pmfW=pmfW,pmfM=pmfM,pmfj=pmfj,gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling, choices=choices)
    }
    jeqfun <- function(theta, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices){
      nl.jacobian(theta, eqfun,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
                  pmfW=pmfW,pmfM=pmfM, pmfj=pmfj,gw=gw,gm=gm,n=n,symmetric=symmetric, sampling=sampling, choices=choices)
    }
    
    ########################
    loglikNULLfun <- function(theta, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices){
      if (symmetric) NumBeta <- dim(Xd)[3] else NumBeta <- dim(Xd)[3]+dim(Zd)[3]
      beta <- rep(0, NumBeta)
      GammaW = theta[1:NumGammaW]
      GammaM = theta[NumGammaW+(1:NumGammaM)]
      delta_w = theta[NumGammaW+NumGammaM+1]
      delta_m = theta[NumGammaW+NumGammaM+2]
      -loglikelihood_CP(beta, GammaW, GammaM, delta_w, delta_m, Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices)
    } 
    eqNULLfun <- function(theta, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices){
      if (symmetric) NumBeta <- dim(Xd)[3] else NumBeta <- dim(Xd)[3]+dim(Zd)[3]
      beta <- rep(0, NumBeta)
      GammaW <- theta[1:NumGammaW]
      GammaM <- theta[NumGammaW+(1:NumGammaM)]
      delta_w = theta[NumGammaW+NumGammaM+1]
      delta_m = theta[NumGammaW+NumGammaM+2]
      equality_constraint_CP(beta, GammaW, GammaM, delta_w, delta_m, Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices)
    }
    gloglikNULLfun <- function(theta, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices){
      nl.grad(theta, loglikNULLfun,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
              pmfW=pmfW,pmfM=pmfM,pmfj=pmfj,gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling, choices=choices)
    }
    jeqNULLfun <- function(theta, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices){
      nl.jacobian(theta, eqNULLfun,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
                  pmfW=pmfW,pmfM=pmfM, pmfj=pmfj,gw=gw,gm=gm,n=n,symmetric=symmetric, sampling=sampling, choices=choices)
    }
    #######################
    
    # for C++ implementation)
  
    # eqfun <- cmpfun(eqfun)
    # loglikfun <- cmpfun(loglikfun)
    # gloglikfun <- cmpfun(gloglikfun)
    # jeqfun <- cmpfun(jeqfun)
    # 
    # eqNULLfun <- cmpfun(eqNULLfun)
    # loglikNULLfun <- cmpfun(loglikNULLfun)
    # gloglikNULLfun <- cmpfun(gloglikNULLfun)
    # jeqNULLfun <- cmpfun(jeqNULLfun)
    
    
    if(control[["algorithm"]]!="solnp"){
        out <- nloptr(x0=theta_0, eval_f=loglikfun, eval_grad_f=gloglikfun,
                      eval_g_eq=eqfun, eval_jac_g_eq=jeqfun,
                      lb=rep(c(-Inf,0,-Inf),c(NumBeta,NumGamma,2)), # 2 for deltas
                      # ub=rep(c(Inf, 1.0e-6),c(NumBeta,NumGamma)),
                      Xd=X,Zd=Z,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
                      pmfW=pmfW, pmfM=pmfM, pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling,
                      choices=choices, opts=control)
        names(out$solution) <- names(theta_0)
        th_hat <- out$solution
        if(!is.null(names(theta_0))){names(th_hat) <- names(theta_0)}
        out$loglik <- -K*out$objective
        out$exitflag <- out$status
        cat("eq values:\n")
        print((eqfun(out$solution, X, Z, NumGammaW, NumGammaM, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices)))
        print(round(th_hat,2))
        
    }else{
        out <- solnp(pars=theta_0, fun=loglikfun, eqfun=eqfun,
                     LB=rep(c(-Inf,0,-Inf),c(NumBeta,NumGamma,2)), # 2 for deltas
                     UB=rep(c(Inf,3000,Inf),c(NumBeta,NumGamma,2)), # 2 for deltas
                     mu=mu,Xd=X,Zd=Z,NumGammaW=NumGammaW,NumGammaM=NumGammaM,
                     pmfj=pmfj,pmfW=pmfW, pmfM=pmfM, choices=choices,
                     control=list(tol=control[["xtol_rel"]],trace=control[["print_level"]],
                                  outer.iter=control[["maxeval"]]))
        names(out$pars) <- names(theta_0)
        th_hat <- out$pars
        out$solution <- out$pars
        #        names(th_hat) <- c(paste("betaW",(1:NumBeta)-1,sep=""), 
        #                           paste("GammaW",(1:NumGammaW)-1,sep=""), paste("GammaM",(1:NumGammaM)-1,sep="") )
        if(!is.null(names(theta_0))){names(th_hat) <- names(theta_0)}
        out$loglik <- -K*out$values[length(out$values)]
        out$exitflag <- out$converged
        cat("eq values:\n")
        print((eqfun(out$solution, X, Z, NumGammaW, NumGammaM, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices)))
        print(round(th_hat,2))
        
    }

    out$eq = eqfun(out$solution, X, Z, NumGammaW, NumGammaM, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling, choices)
    out$coef <- th_hat[1:NumBeta]
    out$rPMFW <- 1/ (1+th_hat[(NumBeta+1):(NumBeta+NumGammaW)])
    out$rPMFM <- 1/ (1+th_hat[(NumBeta+NumGammaW+1):(NumBeta+NumGammaW+NumGammaM)])

    out$loglik <- -K*loglikelihood_CP(th_hat[1:NumBeta],
                                     GammaW=th_hat[NumBeta+(1:NumGammaW)], 
                                     GammaM=th_hat[(NumBeta+NumGammaW)+(1:NumGammaM)],
                                     delta_w = th_hat[NumBeta+NumGammaW+NumGammaM+1],
                                     delta_m = th_hat[NumBeta+NumGammaW+NumGammaM+2],
                                     Xd=X, Zd=Z, pmfW=pmfW, pmfM=pmfM, pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling, choices=choices)

    # null hypothesis
    out.null <- nloptr(x0=theta_0[NumBeta+(1:(NumGammaW+NumGammaM+2))], eval_f=loglikNULLfun, eval_grad_f=gloglikNULLfun,
                  eval_g_eq=eqNULLfun, eval_jac_g_eq=jeqNULLfun,
                  lb=rep(c(0,-Inf),c(NumGamma,2)), # 2 for deltas
                  # ub=rep(c(Inf, 1.0e-6),c(NumBeta,NumGamma)),
                  # Xd=array(0, dim(X)),Zd=array(0, dim(Z)),
                  Xd=X,Zd=Z,
                  NumGammaW=NumGammaW, NumGammaM=NumGammaM,
                  pmfW=pmfW, pmfM=pmfM, pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling,
                  choices=choices, opts=control)
    
    out$loglik.null <- -K*out.null$objective
    #       out$loglik <- out$loglik - out$loglik.null + K*loglik.ref
    #       out$loglik.null <- K*loglik.ref
    out$loglik.null2 <- -K*loglikelihood_CP(rep(0, NumBeta),
                                     GammaW=out.null$solution[1:NumGammaW], 
                                     GammaM=out.null$solution[NumGammaW+(1:NumGammaM)],
                                     delta_w = out.null$solution[NumGamma+1],
                                     delta_m = out.null$solution[NumGamma+2],
                                     Xd=X, Zd=Z, 
                                     pmfW=pmfW, pmfM=pmfM, pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling, choices=choices)
   
    out$null_solution = out.null$solution     
    out$pmfj=pmfj; out$pmfW=pmfW; out$pmfM=pmfM
    out$Xd=X; out$Zd=Z
    
    # chi-squared test
    matching_freq = check_CP_latent(formula, out$solution, mu, Xdata, Zdata, symmetric, choices)
    matching_freq_null = check_CP_latent(formula, c(rep(0, NumBeta), out$null_solution), mu, Xdata, Zdata, symmetric, choices)
    ct = chisq.test(x=matrix(matching_freq$pmfj_est * n,nrow=1)[-(num_Zu+1)*(num_Xu+1)],
               p=matrix(matching_freq_null$pmfj_est,nrow=1)[-(num_Zu+1)*(num_Xu+1)], rescale.p = T)
    out$chisq_stat = ct$statistic
    out$p.value = ct$p.value
                                                
    if(control[["hessian"]]){

      # not centered
      outer_grad = outer_prod_gradient_CP(theta = th_hat, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM,
                                          Xd=X,Zd=Z,pmfW=pmfW, pmfM=pmfM,pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric)
      
      # centered
      outer_grad_2 = outer_prod_gradient_CP_2(theta = th_hat, loglikfun=loglikfun, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM,
                                          Xd=X,Zd=Z,pmfW=pmfW, pmfM=pmfM,pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling)
      
      hess_CP = ave_hessian_CP(theta = th_hat, NumBeta=NumBeta, NumGammaW=NumGammaW, NumGammaM=NumGammaM,
                               Xd=X,Zd=Z,pmfW=pmfW, pmfM=pmfM,pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric)
      
      out$covar = diag(ginv(hess_CP) %*% outer_grad %*% ginv(hess_CP)/n)
      out$covar2 = diag(ginv(hess_CP) %*% outer_grad_2 %*% ginv(hess_CP)/n)
      
      out$outer_grad = outer_grad
      out$outer_grad_2 = outer_grad_2
      out$hess_CP = hess_CP
        
      # browser()
      # dimnames(H) <- list(names(th_hat[1:NumBeta]),names(th_hat[1:NumBeta]))
      # dimnames(Geqfun) <- list(c(names(th_hat)[(NumBeta+1):(NumBeta+NumGamma)],"delta"),names(th_hat[1:NumBeta]))
      # Hi <- try(ginv(-H))
      # if(inherits(Hi,"try-error")){
      #   Hi <- -diag(1/diag(H))
      # }
      # V <- try(Hi - Hi %*% t(Geqfun) %*% ginv(Geqfun %*% Hi %*% t(Geqfun))%*% Geqfun %*% Hi)
      # if(inherits(V,"try-error")){
      #   V <- Hi
      # }
      # if(all(is.na(diag(V)) | abs(diag(V))<1e-8)){
      #   out$covar <- Hi
      # }else{
      #   out$covar <- V
      # }
    }else{
      out$covar <- diag(rep(NA,length(th_hat)))
    }
    
    out$control <- control
    out$df <- K
    
    out$NumBeta <- NumBeta
    out$NumBetaW <- NumBetaW
    out$NumBetaM <- NumBetaM
    out$NumGammaW <- NumGammaW
    out$NumGammaM <- NumGammaM
    out$NumGamma <- NumGamma
    
    class(out) <- "rpm"
    
    return(out)
}
