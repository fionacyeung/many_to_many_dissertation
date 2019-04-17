setwd("C:\\UCLA\\thesis_ideas\\PhD_thesis\\many_to_many_application")

library(knitr)
library(ggplot2)
library(Matrix)

###################### data processing parameters ###############

max_choices = 4 #each panel has 5 wave, got married in wave 2, but the panels are disconnected (can't find the same person in different panels)


library(tictoc)
# library(ggplot2)

source("fitrpm_R_CP.R")
source("loglikelihood_CP.R")
source("equality_constraint_CP.R")
source("rpm.model.matrix.R")
source("choice_probability.R")
source("check_CP.R")
source("asymptotic_var.R")

################## run parameters ############################

# symmetric beta for b1 and b2
symmetric = FALSE

# sampling protocol
sample = "INDIV" # sampling individuals
# sample = "COUPLE" # sampling couples only
# sample = "HOUSEHOLD" # sampling households, which can be a single individual or a couple


################## algorithm parameters #######################

control = list("algorithm"="NLOPT_LD_SLSQP", "symmetric"=symmetric, "sampling_protocol"=sample,
               "xtol_rel"=1.0e-8, "print_level"=2,"maxeval"=1000, 
               "ftol_rel"=1.0e-8,"check_derivatives"=FALSE,"ftol_abs"=1.0e-6,"hessian"=TRUE) 

# starting parameter for estimation
# theta_0 = c(rep(0,numBeta), rep(1,numGamma))
theta_0 = NULL


# ff = ~ b1nodematch("educlevel_t") + b2nodematch("educlevel_t") +
#   b1nodematch("race_t") + b2nodematch("race_t") +
#   b1nodematch("tage_t") + b2nodematch("tage_t")

ff = ~ b1homophily("tage_t") + b2homophily("tage_t") +
  b1greaterthan("tage_t") + b2greaterthan("tage_t") +
  b1smallerthan("tage_t") + b2smallerthan("tage_t")

# Xdata = processed_data$Xdata
# Zdata = processed_data$Zdata
# mu = processed_data$mu

# ff = ~ b1nodematch("educlevel_t") + b2nodematch("educlevel_t")
# out_edu <- fitrpm_R_CP(ff, mu, Xdata, Zdata, theta_0, choices=max_choices, control=control)
# save(out_edu, file="edu_nodematch.RData")
# 
# ff = ~ b1nodematch("race_t") + b2nodematch("race_t")
# out_race <- fitrpm_R_CP(ff, mu, Xdata, Zdata, theta_0, choices=max_choices, control=control)
# save(out_race, file="race_nodematch.RData")
# 
# ff = ~ b1nodematch("tage_t") + b2nodematch("tage_t")
# out_age <- fitrpm_R_CP(ff, mu, Xdata, Zdata, theta_0, choices=max_choices, control=control)
# save(out_age, file="age_nodematch.RData")


Xdata = processed_data_sub$Xdata
Zdata = processed_data_sub$Zdata
mu = processed_data_sub$mu

ff = ~ b1nodematch("educlevel_t") + b2nodematch("educlevel_t")
out_edu_sub <- fitrpm_R_CP(ff, mu, Xdata, Zdata, theta_0, choices=max_choices, control=control)
save(out_edu_sub, file="edu_nodematch.RData")

ff = ~ b1nodematch("race_t") + b2nodematch("race_t")
out_race_sub <- fitrpm_R_CP(ff, mu, Xdata, Zdata, theta_0, choices=max_choices, control=control)
save(out_race_sub, file="race_nodematch.RData")

ff = ~ b1nodematch("tage_t") + b2nodematch("tage_t")
out_age_sub <- fitrpm_R_CP(ff, mu, Xdata, Zdata, theta_0, choices=max_choices, control=control)
save(out_age_sub, file="age_nodematch.RData")



# Compute MLE based on an observed matching
out <- fitrpm_R_CP(ff, mu, Xdata, Zdata, theta_0, choices=max_choices, control=control)
save(out, file="age_homophily_greater_smaller.RData")

print("Coeff:")
print(out$solution)
print("equality:")
print(out$eq)
print("Likelihood ratio:")
print(as.numeric(out$loglik/out$loglik.null))
print("Standard error:")
print(out$covar)
print("Standard error (centered:")
print(out$covar2)

# output table
dff = data.frame(out$solution)
kable(dff)


############ race, edu, age - nodematch ###################
nparam = 1+4+4+4
edu_race_age_coef = out$solution[1:(2*nparam)]
coef_f = edu_race_age_coef[1:nparam]
coef_m = edu_race_age_coef[nparam + 1:nparam]


edu_race_age_df = data.frame(coef = coef_f[-1], 
                             name = c("<HS", "HS", "SomeCollege", "BA+", 
                                      "Hispanic", "Black", "White", "Asian", 
                                      "18-27.9", "28-38.9", "39-48.9", "49-59"),
                             attribute = c(rep("edu", 4), rep("race", 4), rep("age", 4)))

ggplot(edu_race_age_df, aes(x=name, y=coef, fill=attribute)) + geom_bar(stat="identity") +
  scale_x_discrete(limits=edu_race_age_df$name) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("Homophily")


# education by itself
nparam = 5
edu_coef = out_edu$solution[1:(2*nparam)]
edu_coef_f = edu_coef[1:nparam]
edu_coef_m = edu_coef[nparam + 1:nparam]
edu_df = data.frame(coef = c(edu_coef_f[-1], edu_coef_m[-1]), name=rep(c("<HS", "HS", "SomeCollege", "BA+"),2),
                    gender=rep(c("female", "male"), each=nparam-1))
ggplot(edu_df, aes(x=name, y=coef, fill=gender)) + geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(limits=c("<HS", "HS", "SomeCollege", "BA+")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("Homophily")


# race by itself
race_coef = out_race$solution[1:(2*nparam)]
race_coef_f = race_coef[1:nparam]
race_coef_m = race_coef[nparam + 1:nparam]
race_df = data.frame(coef = c(race_coef_f[-1], race_coef_m[-1]), name=rep(c("Hispanic", "Black", "White", "Asian"),2),
                    gender=rep(c("female", "mMale"), each=nparam-1))
ggplot(race_df, aes(x=name, y=coef, fill=gender)) + geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(limits=c("Hispanic", "Black", "White", "Asian")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("Homophily")



# age group by itself
age_coef = out_age$solution[1:(2*nparam)]
age_coef_f = age_coef[1:nparam]
age_coef_m = age_coef[nparam + 1:nparam]
age_df = data.frame(coef = c(age_coef_f[-1], age_coef_m[-1]), name=rep(c("18-27.9", "28-38.9", "39-48.9", "49-59"),2),
                    gender=rep(c("female", "male"), each=nparam-1))
ggplot(age_df, aes(x=name, y=coef, fill=gender)) + geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(limits=c("18-27.9", "28-38.9", "39-48.9", "49-59")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("Homophily")


# # combine edu, race, age plots
# library(gridExtra)
# edu_plot = edu_plot + theme(legend.position = "none") + scale_y_continuous(limits = c(0, 2.5))
# race_plot = race_plot + theme(legend.position = "none") + scale_y_continuous(limits = c(0, 2.5))
# age_plot = age_plot + scale_y_continuous(limits = c(0, 2.5))
# grid.arrange(edu_plot, race_plot, age_plot, ncol=3)

########## age homophily, greater than, smaller than ############

nparam = 4
age_coef = out$solution[1:(2*nparam)]
coef_f = age_coef[1:nparam]
coef_m = age_coef[nparam + 1:nparam]


age_df = data.frame(coef = c(coef_f[-1], coef_m[-1]),
                             name = rep(c("homophily", 
                                      "greater than own", 
                                      "less than own"), 2),
                             gender = c(rep("female", nparam-1), rep("male", nparam-1)))

ggplot(age_df, aes(x=name, y=coef, fill=gender)) + geom_bar(stat="identity", position=position_dodge(width = 0.9)) + # position="dodge") +
  # scale_x_discrete(limits=age_df$name) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("Age Preference")

print(out$covar[1:(nparam*2)])
dff = data.frame(out$solution)
kable(dff)


######### age match type distribution #####################
age_pairs_df = data.frame(pair_id = processed_data$females[,"pair_id"], female_age = processed_data$females[,"tage_t"], 
                          male_age = processed_data$males[,"tage_t"], type = rep(NA, nrow(processed_data$females))) 
uniq_age_pairs = unique(age_pairs_df[,c("female_age", "male_age")])
uniq_age_pairs = uniq_age_pairs[order(uniq_age_pairs$female_age, uniq_age_pairs$male_age),]
uniq_age_pairs$type = 1:nrow(uniq_age_pairs)
  
for(ii in 1:nrow(uniq_age_pairs)) {
  idx = which(age_pairs_df$female_age==uniq_age_pairs$female_age[ii] & age_pairs_df$male_age==uniq_age_pairs$male_age[ii])
  age_pairs_df[idx,"type"] = ii
}
age_pairs_df$female_age = as.factor(age_pairs_df$female_age)
age_pairs_df$male_age = as.factor(age_pairs_df$male_age)
age_pairs_df$type = as.factor(age_pairs_df$type)

ggplot(age_pairs_df, aes(x=type, fill=female_age)) + geom_histogram(position="dodge", stat="count") +
  # scale_x_continuous(breaks=1:16, labels=rep(c("18-27.9", "28-38.9", "39-48.9", "49-59"), 4)) +
  scale_x_discrete(breaks=1:16, labels=rep(c("18-27.9", "28-38.9", "39-48.9", "49-59"), 4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = c(0.85, 0.8)) + 
  scale_fill_discrete(name = "Female Age", labels = c("18-27.9", "28-38.9", "39-48.9", "49-59")) +
  xlab("Male Age") # +
  # theme(axis.text.x = element_text(size=12))

hist_age_pairs = hist(age_pairs_df$type)
min(hist_age_pairs$counts)

############## delete some type 1 marriages (too many). This seems to cause the Hessian to be unstable. #####################################
count_vec = rep(0,nrow(uniq_age_pairs))
for(ii in 1:nrow(uniq_age_pairs)) {
  count_vec[ii] = sum(as.numeric(age_pairs_df$type)==ii)
}
t1_idx = which(age_pairs_df$type==1)
sample_t1 = sample(t1_idx, 300, replace=F)
drop_t1_idx = t1_idx[!(t1_idx %in% sample_t1)]

age_pairs_df_sub = age_pairs_df[-drop_t1_idx,]
sum(age_pairs_df_sub$type==1) # should be 300

age_pairs_df_sub$female_age = as.factor(age_pairs_df_sub$female_age)
age_pairs_df_sub$male_age = as.factor(age_pairs_df_sub$male_age)
age_pairs_df_sub$type = as.factor(age_pairs_df_sub$type)

processed_data_sub = processed_data
processed_data_sub$mu = processed_data_sub$mu[-drop_t1_idx,-drop_t1_idx]
processed_data_sub$Xdata = processed_data_sub$Xdata[-drop_t1_idx,]
processed_data_sub$Zdata = processed_data_sub$Zdata[-drop_t1_idx,]
processed_data_sub$females = processed_data_sub$females[-drop_t1_idx,]
processed_data_sub$males = processed_data_sub$males[-drop_t1_idx,]

age_pairs_df_sub = data.frame(pair_id = processed_data_sub$females[,"pair_id"], female_age = processed_data_sub$females[,"tage_t"], 
                          male_age = processed_data_sub$males[,"tage_t"], type = rep(NA, nrow(processed_data_sub$females))) 

length(which(age_pairs_df_sub$female_age == 1 & age_pairs_df_sub$male_age == 1)) # should be 300

Xdata = processed_data_sub$Xdata
Zdata = processed_data_sub$Zdata
mu = processed_data_sub$mu

ff = ~ b1homophily("tage_t") + b2homophily("tage_t") +
  b1greaterthan("tage_t") + b2greaterthan("tage_t") +
  b1smallerthan("tage_t") + b2smallerthan("tage_t")

# Compute MLE based on an observed matching
out_sub <- fitrpm_R_CP(ff, mu, Xdata, Zdata, theta_0, choices=max_choices, control=control)
save(out_sub, file="age_homophily_greater_smaller_sub.RData")

print("Coeff:")
print(out_sub$solution)
print("equality:")
print(out_sub$eq)
print("Likelihood ratio:")
print(as.numeric(out_sub$loglik/out$loglik.null))
print("Standard error:")
print(sqrt(out_sub$covar))
print("Standard error (centered:")
print(sqrt(out_sub$covar2))

# output table
dff_sub = data.frame(out_sub$solution)
kable(dff_sub)

# ggplot(age_pairs_df_sub, aes(x=type, fill=female_age)) + geom_histogram(position="dodge", stat="count") +
#   # scale_x_continuous(breaks=1:16, labels=rep(c("18-27.9", "28-38.9", "39-48.9", "49-59"), 4)) +
#   scale_x_discrete(breaks=1:16, labels=rep(c("18-27.9", "28-38.9", "39-48.9", "49-59"), 4)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   theme(legend.position = c(0.85, 0.75)) + 
#   scale_fill_discrete(name = "Female Age", labels = c("18-27.9", "28-38.9", "39-48.9", "49-59")) +
#   xlab("Male Age")

# Rename all levels
levels(age_pairs_df_sub$female_age) <- paste0("F:",c("18-27.9", "28-38.9", "39-48.9", "49-59"))
levels(age_pairs_df_sub$male_age) <- paste0("M:",c("18-27.9", "28-38.9", "39-48.9", "49-59"))
ggplot(age_pairs_df_sub, aes(x=male_age)) + geom_histogram(position="dodge", stat="count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ female_age) + xlab("")

########### data distribution for race #######################
race_pairs_df_sub = data.frame(pair_id = processed_data_sub$females[,"pair_id"], female_race = as.factor(processed_data_sub$females[,"race_t"]), 
                              male_race = as.factor(processed_data_sub$males[,"race_t"]), type = rep(NA, nrow(processed_data_sub$females))) 

# Rename all levels
levels(race_pairs_df_sub$female_race) <- paste0("F:",c("Hispanic", "Black", "White", "Asian"))
levels(race_pairs_df_sub$male_race) <- paste0("M:",c("Hispanic", "Black", "White", "Asian"))
ggplot(race_pairs_df_sub, aes(x=male_race)) + geom_histogram(position="dodge", stat="count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ female_race) + xlab("")

########### data distribution for edu #######################
edu_pairs_df_sub = data.frame(pair_id = processed_data_sub$females[,"pair_id"], female_edu = as.factor(processed_data_sub$females[,"educlevel_t"]), 
                               male_edu = as.factor(processed_data_sub$males[,"educlevel_t"]), type = rep(NA, nrow(processed_data_sub$females))) 

# Rename all levels
levels(edu_pairs_df_sub$female_edu) <- paste0("F:",c("<HS", "HS", "SomeCollege", "BA+"))
levels(edu_pairs_df_sub$male_edu) <- paste0("M:",c("<HS", "HS", "SomeCollege", "BA+"))
ggplot(edu_pairs_df_sub, aes(x=male_edu)) + geom_histogram(position="dodge", stat="count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ female_edu) + xlab("")

########## data distribution for age homophily, greater than, and less than ########
age_pairs_df_2 = data.frame(pair_id = processed_data_sub$females[,"pair_id"], female_age = processed_data_sub$females[,"tage_t"], 
                          male_age = processed_data_sub$males[,"tage_t"], female_pair_type = rep(NA, nrow(processed_data_sub$females)),
                          male_pair_type = rep(NA, nrow(processed_data_sub$males))) 

# homophily
age_pairs_df_2$female_pair_type[age_pairs_df_2$female_age == age_pairs_df_2$male_age] = "homophily"
age_pairs_df_2$male_pair_type[age_pairs_df_2$female_age == age_pairs_df_2$male_age] = "homophily"

# greater than own
age_pairs_df_2$female_pair_type[age_pairs_df_2$female_age < age_pairs_df_2$male_age] = "greater than own"
age_pairs_df_2$male_pair_type[age_pairs_df_2$male_age < age_pairs_df_2$female_age] = "greater than own"

# less than own
age_pairs_df_2$female_pair_type[age_pairs_df_2$female_age > age_pairs_df_2$male_age] = "less than own"
age_pairs_df_2$male_pair_type[age_pairs_df_2$male_age > age_pairs_df_2$female_age] = "less than own"

ggplot(age_pairs_df_2, aes(x=female_pair_type)) + geom_histogram(position="dodge", stat="count") +
  xlab("Female Age Preference")

ggplot(age_pairs_df_2, aes(x=male_pair_type)) + geom_histogram(position="dodge", stat="count") +
  xlab("Male Age Preference")

age_pairs_df_3 = data.frame(pair_type = c(age_pairs_df_2$female_pair_type, age_pairs_df_2$male_pair_type),
                            gender = rep(c("female", "male"), each=nrow(age_pairs_df_2)))

ggplot(age_pairs_df_3, aes(x=pair_type)) + geom_histogram(position="dodge", stat="count") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(. ~ gender) + xlab("")

##################################

load("out_sub_age_homophily_greater_smaller.Rdata")
print(signif(sqrt(out_sub$covar2),digits=4))
dff_sub = signif(out_sub$solution,digits=4)
kable(dff_sub)

###############################
load("processed_data_sub_age_homophily_greater_smaller.Rdata")
ggplot(data.frame(length=factor(processed_data_sub$mu@x)), aes(length))+geom_histogram(stat="count")+xlab("length of union")
length_unions = hist(processed_data_sub$mu@x)$counts

