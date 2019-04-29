# setwd("C:\\UCLA\\research_projects\\many_to_many_dissertation")
# source("unpartnered_partnered_3.R")

set.seed(1234)

library(knitr)
library(ggplot2)
# library(entropy)
library(latex2exp)

###################### files and libraries for the data ######################

#setwd("C://UCLA//thesis_ideas//homeless_poverty//PAA_annual_meeting_2018//poster//unpartnered_partnered")
data_file = "C:\\UCLA\\thesis_ideas\\homeless_poverty\\Unpartnered&NewPartnered_Final\\Unpartnered&NewPartnered_Final.dta"

library(Matrix)
library(readstata13)

data = read.dta13(data_file)


###################### data processing parameters ###############

complete_cases_only = TRUE

# only use marriages occured in wave 2 of any given panel
# hh=hist(paired$wavetime) # 4198 occurrences in wave 2
wave_num = 2 # only consider marriages that start at this wave
max_choices = 4 #each panel has 5 wave, got married in wave 2, but the panels are disconnected (can't find the same person in different panels)

###################### sanity check for the data #####################

idx = which(!is.na(data$allrespartid_w))
paired = data[idx,c("pid", "spanel_all", "wpfinwgt_t", "tage_t", "female_t", "race_t", "educlevel_t", "allrespartid_w", "newrelunpar" )]

# sanity check
oddidx = seq(1,nrow(paired),by=2)
odds=paired[oddidx,]
evens = paired[-oddidx,]
consq = as.numeric(row.names(odds))-as.numeric(row.names(evens))
all(consq==-1) # this should be true since couples are on consecutive rows
nrow(paired)%%2 == 0 # pairs are paired
unique(data$spanel_all) # unique panels
all(paired$newrelunpar==1) # all pairs are new (newrelupar = 1)

################## select data to include #######################

process_data = function(data, complete_cases_only, wave_num, max_choices) {

  # separate the couples from the singles
  idx = which(!is.na(data$allrespartid_w))
  paired = data[idx,c("pid", "epppnum", "spanel_all", "wavetime", "wpfinwgt_t", "tage_t", "female_t", "race_t", "educlevel_t", "allrespartid_w", "newrelunpar" )]
  paired$pair_id = rep(1:(nrow(paired)/2), each=2)
  single = data[-idx,c("pid", "epppnum", "spanel_all", "wavetime", "wpfinwgt_t", "tage_t", "female_t", "race_t", "educlevel_t", "allrespartid_w", "newrelunpar" )]
  single$pair_id = rep(0,nrow(single))
  
  # order the partners at the same index
  paired_females = paired[paired$female_t == "Female",] # their spouses have the same indices in paired_males
  single_females = single[single$female_t == "Female",]
  paired_males = paired[paired$female_t == "Male",]
  single_males = single[single$female_t == "Male",]
  # sanity check
  nrow(paired_females) == nrow(paired_males)
  
  paired2 = paired[paired$wavetime==wave_num,] # 4198 marriages
  single2 = single[single$pid %in% paired2$pid & single$wavetime > wave_num,] # 444 divorces (~30% divorce rate)
  
  # get the people who remarried
  remarried_idx = duplicated(paired$pid)
  remarried = paired[remarried_idx,]
  remarried_from_w2 = remarried[remarried$pid %in% paired2$pid,] # divorced from wave 2 and married again later in the panel
  num_times_remarried2 = as.data.frame(table(remarried_from_w2$pid))
  # any(duplicated(remarried_from_w2$pid)) # no one remarried 3 times
  
  # get the length of marriages from wave 2
  num_times_single2 = as.data.frame(table(single2$pid))
  temp = merge(num_times_single2, num_times_remarried2, by="Var1", all=T) # outer join the remarry and single reports on the same person
  temp$Freq.x[is.na(temp$Freq.x)]=0
  temp$Freq.y[is.na(temp$Freq.y)]=0
  temp$single = temp$Freq.x+temp$Freq.y
  # max_choices - temp$total # number of waves the marriage from wave 2 lasted (but this can be wrong if the second marriage lasts till wave 5 -- minority case)
  
  all(temp$Var1 %in% paired2$pid ) # sanity check
  temp = temp[,c("Var1", "single")]
  paired2$paired = max_choices
  paired2 = merge(paired2, temp, by.x="pid", by.y="Var1", all=T) # outer joint
  paired2$single[is.na(paired2$single)] = 0
  paired2$paired = paired2$paired - paired2$single
    
  # order by the pair_id again (was scrambled by merge)
  paired2 = paired2[order(paired2$pair_id, paired2$paired, -paired2$single),]
  # fix the problem where secondary correspondent not reported as single after divorce
  primary = paired2$paired[seq(1, nrow(paired2) , by=2)] # primary correspondent's paired value
  paired2$paired = rep(primary, each=2)
  primary = paired2$single[seq(1, nrow(paired2) , by=2)] # primary correspondent's single value
  paired2$single = rep(primary, each=2)
  
  # split into females and males
  females = paired2[paired2$female_t == "Female",] # their spouses have the same indices in paired_males
  males = paired2[paired2$female_t == "Male",]
  
  if (complete_cases_only) {
    # get rid of the incomplete cases for pairs
    incompl_idx = which(!complete.cases(females))
    females = females[-incompl_idx,]
    males = males[-incompl_idx,]
    incompl_idx = which(!complete.cases(males))
    females = females[-incompl_idx,]
    males = males[-incompl_idx,]

    # sanity check
    nrow(females) == nrow(males) # 1953 females and 1953 males left
  }
  
  # get the weight for the pair (if both partners are from the same wave, just use the one with lower epppnum)
  pair_w_idx = apply(cbind(females$epppnum, males$epppnum), 1, which.min)
  pair_w = cbind(females$wpfinwgt_t, males$wpfinwgt_t)
  pair_w = pair_w[cbind(1:length(pair_w_idx),pair_w_idx)]
  pair_min_epppnum = apply(cbind(females$epppnum, males$epppnum), 1, min)
  # sanity check
  all(pair_min_epppnum < 200)
  
  Xdata = females
  Zdata = males
  
  # sanity checks
  all(Xdata$female_t == "Female")
  all(Zdata$female_t == "Male")
  
  
  # only keep columns that are attributes
  Xdata = Xdata[,c("tage_t", "race_t", "educlevel_t")]
  Zdata = Zdata[,c("tage_t", "race_t", "educlevel_t")]
  
  
  # create the adjacency matrix
  n_pairs = nrow(Xdata)

  # observed matching (sparse matrix)
  mu = Diagonal(n_pairs) * females$paired
  
  Xdata = as.matrix(Xdata)
  Zdata = as.matrix(Zdata)
  
  return(list(Xdata=Xdata, Zdata=Zdata, mu=mu, 
              females=females, males=males, 
              X_w=females$wpfinwgt_t, Z_w=males$wpfinwgt_t , pair_w=pair_w))
}




######################## pre-process data ##############################

# copy co-habitant's characteristics into the spouse's characteristics 
# (we don't distinguish cohabitation and marriage for now)
cohabidx = which(data$marrcohabt_rev==2)
data$s_wpfinwgt_t[cohabidx] = data$p_wpfinwgt_t[cohabidx]
data$s_tage_t[cohabidx] = data$p_tage_t[cohabidx]
data$s_female_t[cohabidx] = data$p_female_t[cohabidx]
data$s_race_t[cohabidx] = data$p_race_t[cohabidx]
data$s_educlevel_t[cohabidx] = data$p_educlevel_t[cohabidx]

# use numeric values for factors
data$educlevel_t = as.character(data$educlevel_t)
data$educlevel_t[data$educlevel_t == "<HS"] = 1
data$educlevel_t[data$educlevel_t == "HS"] = 2
data$educlevel_t[data$educlevel_t == "SomeCollege"] = 3
data$educlevel_t[data$educlevel_t == "BA+"] = 4
data$educlevel_t = as.numeric(data$educlevel_t)
data$race_t = as.character(data$race_t)
data$race_t[data$race_t == "Hispanic"] = 1
data$race_t[data$race_t == "Black"] = 2
data$race_t[data$race_t == "White"] = 3
data$race_t[data$race_t == "Asian"] = 4
data$race_t = as.numeric(data$race_t)

race_label = c("Hispanic", "Black", "White", "Asian")
age_label = levels(cut(data$tage_t, breaks = 4))
data$tage_t <- cut(data$tage_t, breaks = 4, labels=FALSE)


###################### more data processing ###################

processed_data = process_data(data, complete_cases_only, wave_num, max_choices)

########## delete some type 1 marriages (too many). This seems to cause the Hessian to be unstable. ########

# age_pairs_df = data.frame(pair_id = processed_data$females[,"pair_id"], female_age = processed_data$females[,"tage_t"],
#                           male_age = processed_data$males[,"tage_t"], type = rep(NA, nrow(processed_data$females)))
# uniq_age_pairs = unique(age_pairs_df[,c("female_age", "male_age")])
# uniq_age_pairs = uniq_age_pairs[order(uniq_age_pairs$female_age, uniq_age_pairs$male_age),]
# uniq_age_pairs$type = 1:nrow(uniq_age_pairs)
# 
# for(ii in 1:nrow(uniq_age_pairs)) {
#   idx = which(age_pairs_df$female_age==uniq_age_pairs$female_age[ii] & age_pairs_df$male_age==uniq_age_pairs$male_age[ii])
#   age_pairs_df[idx,"type"] = ii
# }
# age_pairs_df$female_age = as.factor(age_pairs_df$female_age)
# age_pairs_df$male_age = as.factor(age_pairs_df$male_age)
# age_pairs_df$type = as.factor(age_pairs_df$type)
# 
# count_vec = rep(0,nrow(uniq_age_pairs))
# for(ii in 1:nrow(uniq_age_pairs)) {
#   count_vec[ii] = sum(as.numeric(age_pairs_df$type)==ii)
# }
# t1_idx = which(age_pairs_df$type==1)
# sample_t1 = sample(t1_idx, 300, replace=F)
# drop_t1_idx = t1_idx[!(t1_idx %in% sample_t1)]
# 
# age_pairs_df_sub = age_pairs_df[-drop_t1_idx,]
# sum(age_pairs_df_sub$type==1) # should be 300
# 
# age_pairs_df_sub$female_age = as.factor(age_pairs_df_sub$female_age)
# age_pairs_df_sub$male_age = as.factor(age_pairs_df_sub$male_age)
# age_pairs_df_sub$type = as.factor(age_pairs_df_sub$type)
# 
# processed_data_sub = processed_data
# processed_data_sub$mu = processed_data_sub$mu[-drop_t1_idx,-drop_t1_idx]
# processed_data_sub$Xdata = processed_data_sub$Xdata[-drop_t1_idx,]
# processed_data_sub$Zdata = processed_data_sub$Zdata[-drop_t1_idx,]
# processed_data_sub$females = processed_data_sub$females[-drop_t1_idx,]
# processed_data_sub$males = processed_data_sub$males[-drop_t1_idx,]
# processed_data_sub$X_w = processed_data_sub$X_w[-drop_t1_idx]
# processed_data_sub$Z_w = processed_data_sub$Z_w[-drop_t1_idx]
# processed_data_sub$pair_w = processed_data_sub$pair_w[-drop_t1_idx]
# 
# age_pairs_df_sub = data.frame(pair_id = processed_data_sub$females[,"pair_id"], female_age = processed_data_sub$females[,"tage_t"],
#                               male_age = processed_data_sub$males[,"tage_t"], type = rep(NA, nrow(processed_data_sub$females)))
# 
# length(which(age_pairs_df_sub$female_age == 1 & age_pairs_df_sub$male_age == 1)) # should be 300
# 
# Xdata = processed_data_sub$Xdata
# Zdata = processed_data_sub$Zdata
# mu = processed_data_sub$mu
# 
# save(processed_data_sub, file="processed_data_sub_weight.RData")

# load("processed_data_sub.Rdata")
load("processed_data_sub_weight.RData")

processed_data = processed_data_sub

print(paste0("# females = ", nrow(processed_data$females), ", # males = ", nrow(processed_data$males)))

################# source the algo files and load the libraries ################

library(tictoc)
# library(ggplot2)

source("fitrpm_R_CP.R")
source("loglikelihood_CP.R")
source("equality_constraint_CP.R")
source("rpm.model.matrix.R")
source("choice_probability.R")
source("check_CP.R")
source("asymptotic_var.R")

################## source some files for plotting data ################

# source("edu_hypo_hyper_same.R")

################## set forumula ##########################

# specify the formula for utilities
# example: ff = ~ b1cov("f1") + b2cov("f1") + b1absdiff("f1",1) + b2absdiff("f1",1)


# used in dissertation
ff = ~ b1nodematch("educlevel_t") + b2nodematch("educlevel_t") # used in dissertation

# ff = ~ b1nodematch("race_t") + b2nodematch("race_t")

# ff = ~ b1nodematch("tage_t") + b2nodematch("tage_t")

# ff = ~ b1homophily("tage_t") + b2homophily("tage_t") +
#   b1greaterthan("tage_t") + b2greaterthan("tage_t") +
#   b1smallerthan("tage_t") + b2smallerthan("tage_t")



# takes > 1 min
# ff = ~ b1homophily("educlevel_t") + b2homophily("educlevel_t") + b1homophily("race_t") + b2homophily("race_t")
# ff = ~ b1nodematch("educlevel_t") + b2nodematch("educlevel_t") + b1homophily("race_t") + b2homophily("race_t") # takes a a couple minutes



# ff = ~ b1nodematch("educlevel_t") + b2nodematch("educlevel_t") +
#   b1nodematch("tage_t") + b2nodematch("tage_t")

# takes too long (don't run this!)
# ff = ~ b1nodematch("educlevel_t") + b2nodematch("educlevel_t") +
#   b1nodematch("race_t") + b2nodematch("race_t") +
#   b1nodematch("tage_t") + b2nodematch("tage_t")

# ff = ~ b1nodematch("educlevel_t") + b2nodematch("educlevel_t") +
#   b1nodematch("race_t") + b2nodematch("race_t") 

# interactions
# ff = ~ b1homophily("educlevel_t","race_t") + b2homophily("educlevel_t","race_t")

################ parallel processing ################

# # parallel
# cl<-makeCluster(4) #change to your number of CPU cores
# # cl <- makeCluster(16, type="MPI")
# registerDoSNOW(cl)

reme = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

################## run parameters ############################

# symmetric beta for b1 and b2
symmetric = FALSE

# sampling protocol
sample = "INDIV" # sampling individuals
# sample = "COUPLE" # sampling couples only
# sample = "HOUSEHOLD" # sampling households, which can be a single individual or a couple


################## algorithm parameters #######################

control = list("algorithm"="NLOPT_LD_SLSQP", "symmetric"=symmetric, "sampling_protocol"=sample,
               "xtol_rel"=1.0e-8, "print_level"=0,"maxeval"=1000, 
               "ftol_rel"=1.0e-8,"check_derivatives"=FALSE,"ftol_abs"=1.0e-6,"hessian"=TRUE) 

# starting parameter for estimation
# theta_0 = c(rep(0,numBeta), rep(1,numGamma))
theta_0 = NULL

################################################################
################# run algorithm ################################



tic.clearlog()
tic("start")

############## bootstrap ###################

library(foreach)
library(doSNOW)

# parallel
cl<-makeCluster(4) #change to your number of CPU cores
# cl <- makeCluster(16, type="MPI")
registerDoSNOW(cl)

# out$solution, out$eq, w_ave_choices, m_ave_choices, out$null_solution, out$chisq_stat, out$p.value, out$covar2
# bootstrap_result <- matrix(0,ncol=(numBeta+2+numGamma)+(numGamma+1)+2+(numGamma+2)+2+(numBeta+numGamma+2),nrow=B) 

B = 1000
# B = 3

bootstrap_result <-
 foreach (b=1:B, .combine = 'rbind', .packages=c('nloptr','abind', 'Matrix', 'numDeriv', 'MASS', 'questionr')) %dopar% {
# for (b in 1:B) {    
   
    keep_idx = sample(1:nrow(processed_data$mu), nrow(processed_data$mu), replace=T)
    
    ############# fit data #####################
    
    Xdata = processed_data$Xdata[keep_idx,]
    Zdata = processed_data$Zdata[keep_idx,]
    X_w = processed_data$X_w[keep_idx]
    Z_w = processed_data$Z_w[keep_idx]
    pair_w = processed_data$pair_w[keep_idx]
    mu = Diagonal(nrow(processed_data$mu)) * diag(processed_data$mu)[keep_idx]
    
    # Compute MLE based on an observed matching
    out <- fitrpm_R_CP(ff, mu, Xdata, Zdata, X_w, Z_w, pair_w, theta_0, choices=max_choices, control=control)
    # save(out, file="edu_race_age_nodematch.RData")
    # c(out$solution, out$eq, w_ave_choices, m_ave_choices, out$null_solution, out$chisq_stat, out$p.value, out$covar2)
    c(out$solution, out$eq, out$null_solution, out$chisq_stat, out$p.value, out$covar2)
    
  } # end of bootstrap loop
toc(log=TRUE, quiet=FALSE)

stopCluster(cl)

# original data set
Xdata = processed_data_sub$Xdata
Zdata = processed_data_sub$Zdata
mu = processed_data_sub$mu
X_w = processed_data_sub$X_w
Z_w  = processed_data_sub$Z_w
pair_w = processed_data_sub$pair_w

out <- fitrpm_R_CP(ff, mu, Xdata, Zdata, X_w, Z_w, pair_w, theta_0, choices=max_choices, control=control)

print("coeff:")
print(out$solution)
print("equality:")
print(out$eq)
print("Likelihood ratio:")
print(as.numeric(out$loglik/out$loglik.null))
print("Standard error:")
print(sqrt(out$covar))
print("Standard error (centered:")
print(sqrt(out$covar2))

# output table
dff = data.frame(out$solution)
print(kable(dff))

# compute variance and standard error of the estimates
diff_est = bootstrap_result[,1:length(out$solution)] - matrix(rep(out$solution, times=B), nrow=B, byrow=T)
asympt_var = matrix(0, nrow=length(out$solution), ncol=length(out$solution))
for (b in 1:B) {
  asympt_var = asympt_var + outer(diff_est[b,], diff_est[b,])
}
asympt_var = asympt_var/B
se = sqrt(diag(asympt_var))

save.image(file=paste0("SE_B", B, ".RData"))
print("bootstrap standard error:")
print(se)

# compare estimated joint probabilities with truth
pmfj = check_CP_latent(ff, out$solution, mu, Xdata, Zdata, symmetric, max_choices)
print("estimated joint probabilities")
print(pmfj$pmfj_est)
print("observed joint probabilities")
print(pmfj$pmfj_obs)


log.txt <- tic.log(format = TRUE)
log.lst <- tic.log(format = FALSE)
tic.clearlog()
timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
writeLines(unlist(log.txt))


############### plot the data ####################

library(ggplot2)

# plot by education
joint_edu_df = data.frame(self_type = processed_data$Xdata[,"educlevel_t"], spouse_type = processed_data$Zdata[,"educlevel_t"])
# joint_edu_df = joint_edu_df[rep(row.names(joint_edu_df), processed_data$females[,"paired"]),]
joint_edu_df$self_type = as.factor(joint_edu_df$self_type)
# joint_edu_df$spouse_type = as.factor(joint_edu_df$spouse_type)

ggplot(joint_edu_df, aes(x=spouse_type, fill=self_type)) +  geom_histogram(binwidth=1, alpha=.5, position="dodge") + 
  scale_x_discrete(name ="Spouse Education", limits=c("<HS","HS","SomeCollege","BA+")) + 
  ggtitle("Education Homophily") + scale_fill_discrete(name="Self Education",
                                                        #breaks=c("1", "2", "3", "4"),
                                                        labels=c("<HS","HS","SomeCollege","BA+"))


# plot by age group
joint_age_df = data.frame(self_type = processed_data$Xdata[,"tage_t"], spouse_type = processed_data$Zdata[,"tage_t"])
# joint_age_df = joint_age_df[rep(row.names(joint_age_df), processed_data$females[,"paired"]),]
joint_age_df$self_type = as.factor(joint_age_df$self_type)
# joint_age_df$spouse_type = as.factor(joint_age_df$spouse_type)

ggplot(joint_age_df, aes(x=spouse_type, fill=self_type)) +  geom_histogram(binwidth=1, alpha=.5, position="dodge") + 
  scale_x_discrete(name ="Spouse Age Group", limits=age_label) + 
  ggtitle("Age Homophily") + scale_fill_discrete(name="Self Age",
                                                        #breaks=c("1", "2", "3", "4"),
                                                        labels=age_label)


# plot by race
joint_race_df = data.frame(self_type = processed_data$Xdata[,"race_t"], spouse_type = processed_data$Zdata[,"race_t"])
# joint_race_df = joint_race_df[rep(row.names(joint_race_df), processed_data$females[,"paired"]),]
joint_race_df$self_type = as.factor(joint_race_df$self_type)
# joint_race_df$spouse_type = as.factor(joint_race_df$spouse_type)

ggplot(joint_race_df, aes(x=spouse_type, fill=self_type)) +  geom_histogram(binwidth=1, alpha=.5, position="dodge") + 
  scale_x_discrete(name ="Spouse Race", limits=race_label) + 
  ggtitle("Race Homophily") + scale_fill_discrete(name="Self Race",
                                                 #breaks=c("1", "2", "3", "4"),
                                                 labels=race_label)


