# setwd("C:\\UCLA\\research_projects\\many_to_many_dissertation")
# source("variances_of_observed_and_unobserved_factors.R")

library(tictoc)
library(ggplot2)
library(Matrix)
library(foreach)
library(doSNOW)
library(abind)
library(latex2exp)
library(knitr)

source("rpm.model.matrix.R")

set.seed(1234)
# set.seed(712434)

# parallel
cl<-makeCluster(4) #change to your number of CPU cores
# cl <- makeCluster(16, type="MPI")
registerDoSNOW(cl)

reme = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

################## create simulated data ##########################

############## temp to match one-to-one #####################

# specify the formula for utilities
ff = ~ b1cov("f1") + b2cov("f1") + b1absdiff("f1",1) + b2absdiff("f1",1)
# ff = ~ b1cov("f1") + b2cov("f1") + b1cov("f2") + b2cov("f2")

#############################################################

# parse the formula
model.terms <- rownames(attr(terms.formula(ff), "factors"))
temp <- strsplit(model.terms, "[(]")
model.terms.names <- unlist(lapply(temp, `[[`, 1))
temp2 <- gsub(")", "", gsub("\"", "", unlist(temp)))
temp2 <- gsub("[ \t\n\r\f\v]", "", temp2)
model.terms.coef.names <- temp2[-which(temp2 %in% model.terms.names)]
model.terms.coef.names <- strsplit(model.terms.coef.names, ",")


# true parameter for utility generation
# beta_wt <- c(1,0,1)
# beta_mt <- c(1,0,1)
beta_wt = c(0.5,0.5,1)
beta_mt = c(0.5,0.5,1)
beta_w <- matrix(beta_wt,ncol=1)
beta_m <- matrix(beta_mt,ncol=1)

############# simulation parameters and result ####################

N <- c(20,50,100,200,500,1000,2000)
# N = c(20, 50)
B = 600

########## modified to work with many-to-many ###################
# choices = 1:6
############################################################


# portions of men and women
gw = log(0.5*2) # 50% women: exp(gw)+exp(gm) has to be 2
gm = log(0.5*2) # 50% men


var_sim <-
  # foreach(cc=choices, .combine = 'rbind') %:%
  foreach(n=N, .combine = 'rbind') %:%
    foreach (b=1:B, .combine = 'rbind', .packages=c('nloptr','abind', 'Matrix')) %dopar% {
      
      print(paste0("N=",n, " b=",b))
      
      ######### temp to match one-to-one ###############
      nw = round(n*exp(gw))
      nm = round(n*exp(gm))
      
      # create X, Z data
      X = matrix(rep(0:1,times=nw/2),ncol=1)
      colnames(X) = "f1"
      Z = matrix(rep(0:1,times=nm/2),ncol=1)
      colnames(Z) = "f1"
      
      
      ##################################################
      
      
      # create utility matrices
      # U_star = matrix(0, nrow=nw, ncol = nm)
      U_star = matrix(0, nrow=dim(X)[1], ncol = dim(Z)[1])
      # V_star = matrix(0, nrow=nm, ncol = nw)
      V_star = matrix(0, nrow=dim(Z)[1], ncol = dim(X)[1])
      modelmat <- rpm.model.matrix(model.terms.names, model.terms.coef.names, X, Z)
      for (ii in 1:dim(modelmat$X)[3]) {
        U_star = U_star + modelmat$X[,,ii] * beta_w[ii]
      }
      for (ii in 1:dim(modelmat$Z)[3]) {
        V_star = V_star + modelmat$Z[,,ii] * beta_m[ii]
      }
      # Menzel's code did t(V_star) -- each column is a man (same as in U_star)
      V_star = t(V_star)
      
      
      var_U_star = var(matrix(U_star, ncol=1))
      var_V_star = var(matrix(V_star, ncol=1))
      
      
      # eta  <- -log(-log(matrix(runif(nw * nm), nw)))
      # zeta <- -log(-log(matrix(runif(nw * nm), nw)))
      eta  <- -log(-log(matrix(runif(dim(X)[1] * dim(Z)[1]), dim(X)[1])))
      zeta <- -log(-log(matrix(runif(dim(X)[1] * dim(Z)[1]), dim(X)[1])))
      
      var_eta_paired = var(matrix(eta, ncol=1))
      var_zeta_paired = var(matrix(zeta, ncol=1))
      
      mean_eta_paired = mean(matrix(eta, ncol=1))
      mean_zeta_paired = mean(matrix(zeta, ncol=1))
      
      # adjust for outside option (remain single)
      J <- round(sqrt(n))
      # eta0 <- reme(matrix(apply(-log(-log(matrix(runif(J * nw), J))), 2, max),ncol=1),1,nm)
      # zeta0 <- reme(matrix(apply(-log(-log(matrix(runif(J * nm), J))),2,max),nrow=1), nw, 1)
      eta0 <- reme(matrix(apply(-log(-log(matrix(runif(J * dim(X)[1]), J))), 2, max),ncol=1),1,dim(Z)[1])
      zeta0 <- reme(matrix(apply(-log(-log(matrix(runif(J * dim(Z)[1]), J))),2,max),nrow=1), dim(X)[1], 1)
      eta <- eta - eta0
      zeta <- zeta - zeta0
      
      var_eta0 = var(matrix(eta0, ncol=1))
      var_zeta0 = var(matrix(zeta0, ncol=1)) 
      
      mean_eta0 = mean(matrix(eta0, ncol=1))
      mean_zeta0 = mean(matrix(zeta0, ncol=1)) 
      
      var_eta = var(matrix(eta, ncol=1))
      var_zeta = var(matrix(zeta, ncol=1)) 
      
      mean_eta = mean(matrix(eta, ncol=1))
      mean_zeta = mean(matrix(zeta, ncol=1)) 
      
      U <- U_star + eta
      V <- V_star + zeta
      
      var_U = var(matrix(U, ncol=1))
      var_V = var(matrix(V, ncol=1))
      
      # c(n, cc, b, var_U_star, var_V_star, var_eta_paired, var_zeta_paired, var_eta0, var_zeta0, var_eta, var_zeta, var_U, var_V)
      c(n, b, var_U_star, var_V_star, var_eta_paired, var_zeta_paired, var_eta0, var_zeta0, var_eta, var_zeta, var_U, var_V,
        mean_eta_paired, mean_zeta_paired, mean_eta0, mean_zeta0, mean_eta, mean_zeta)
      
    } # end of B loop


stopCluster(cl)

# colnames(var_sim) = c("N", "choice", "b", "U_star", "V_star", "eta_paired", "zeta_paired", "eta_unpaired", "zeta_unpaired", 
#                       "eta", "zeta", "U", "V", "m_eta_paired", "m_zeta_paired", "m_eta_unpaired", "m_zeta_unpaired", 
#                       "m_eta", "m_zeta")
colnames(var_sim) = c("N", "b", "U_star", "V_star", "eta_paired", "zeta_paired", "eta_unpaired", "zeta_unpaired", 
                      "eta", "zeta", "U", "V", "m_eta_paired", "m_zeta_paired", "m_eta_unpaired", "m_zeta_unpaired", 
                      "m_eta", "m_zeta")

df = as.data.frame(var_sim)
df$N = as.factor(df$N)
df$choice = as.factor(df$choice)
df$b = as.factor(df$b)

df$var_error = df$U / df$eta - 1

# appender <- function(string) {
#   # if (string[1]=="consumer" || string[1]=="supplier")
#   #   string
#   # else
#     TeX(paste0("$\\tau = $", string)) 
# }
# 
# appender2 <- function(string) {
#   # if (string[1]=="consumer" || string[1]=="supplier")
#   #   string
#   # else
#   TeX(paste0("$\\N = $", string)) 
# }


# print(ggplot(df, aes(U, colour=N)) + geom_density() + facet_grid(. ~ choice, labeller = as_labeller(appender, default = label_parsed))) 
# print(ggplot(df, aes(eta, colour=N)) + geom_density() + facet_grid(. ~ choice, labeller = as_labeller(appender, default = label_parsed)))

# probability density function for the ratio of variances of observed to unobserved factors
print(ggplot(df, aes(var_error, colour=N)) + geom_density() + xlab("ratio of variances") + labs(color = "n"))

mean_var = tapply(df$var_error, df$N, mean)
print("mean of variance ratio:")
print(mean_var)

var_var = tapply(df$var_error, df$N, var)
print("variance of variance ratio:")
print(var_var)

mean_eta_p = tapply(df$m_eta_paired, df$N, mean)
mean_eta_up = tapply(df$m_eta_unpaired, df$N, mean)
mean_eta = tapply(df$m_eta, df$N, mean)

df2 = data.frame(mean_ratio = mean_var, var_ratio = var_var, 
                 mean_eta_p=mean_eta_p, mean_eta_up=mean_eta_up, mean_eta=mean_eta)
kable(df2)

save.image("variances.RData")


