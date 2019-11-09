# setwd("C:\\UCLA\\research_projects\\many_to_many_dissertation")
# source("survival_analysis.R")

# see tutorials at https://www.datacamp.com/community/tutorials/survival-analysis-R
# and http://www.sthda.com/english/wiki/cox-proportional-hazards-model

library(survival)
library(survminer)
library(Matrix)

######### tutorial example ###################

# data(ovarian)
# 
# ovarian$rx <- factor(ovarian$rx, 
#                      levels = c("1", "2"), 
#                      labels = c("A", "B"))
# ovarian$resid.ds <- factor(ovarian$resid.ds, 
#                            levels = c("1", "2"), 
#                            labels = c("no", "yes"))
# ovarian$ecog.ps <- factor(ovarian$ecog.ps, 
#                           levels = c("1", "2"), 
#                           labels = c("good", "bad"))
# 
# ovarian$age_group = "young"
# ovarian$age_group[ovarian$age >= 50] = "old"
# ovarian$age_group <- factor(ovarian$age_group)
# 
# # Fit survival data using the Kaplan-Meier method
# surv_object <- Surv(time = ovarian$futime, event = ovarian$fustat) # futime = time of death or censored, fustat = censored or not  
# # surv_object
# 
# fit1 <- survfit(surv_object ~ rx, data = ovarian)
# summary(fit1)
# ggsurvplot(fit1, data = ovarian, pval = TRUE)
# 
# fit2 <- survfit(surv_object ~ resid.ds, data = ovarian)
# ggsurvplot(fit2, data = ovarian, pval = TRUE)
# 
# 
# # Fit a Cox proportional hazards model
# fit.coxph <- coxph(surv_object ~ rx + resid.ds + age_group + ecog.ps, 
#                    data = ovarian)
# ggforest(fit.coxph, data = ovarian)


########################################################

# set.seed(1234)

load("processed_data_sub.Rdata")

# use one variable to indicate age homphily, greater than own, or less than own
age_disparity_f = matrix(0,nrow=nrow(processed_data_sub$Xdata),ncol=1)
idx_greater = as.numeric(which((processed_data_sub$Xdata[,"tage_t"] ==1 & 
                                  processed_data_sub$Zdata[,"tage_t"] > 1) |
                                 (processed_data_sub$Xdata[,"tage_t"] ==2 & 
                                    processed_data_sub$Zdata[,"tage_t"] > 2) | 
                                 (processed_data_sub$Xdata[,"tage_t"] ==3 & 
                                    processed_data_sub$Zdata[,"tage_t"] > 3)))
age_disparity_f[idx_greater] = 1

idx_less = as.numeric(which((processed_data_sub$Xdata[,"tage_t"] ==4 & 
                               processed_data_sub$Zdata[,"tage_t"] < 4) |
                              (processed_data_sub$Xdata[,"tage_t"] ==3 & 
                                 processed_data_sub$Zdata[,"tage_t"] < 3) | 
                              (processed_data_sub$Xdata[,"tage_t"] ==2 & 
                                 processed_data_sub$Zdata[,"tage_t"] < 2)))
age_disparity_f[idx_less] = 2


# nodematch edu
homophily_edu_f = matrix(0,nrow=nrow(processed_data_sub$Xdata),ncol=1)
homophily_edu_f[as.numeric(which(processed_data_sub$Xdata[,"educlevel_t"]==1 & 
         processed_data_sub$Zdata[,"educlevel_t"]==1))] = 1
homophily_edu_f[as.numeric(which(processed_data_sub$Xdata[,"educlevel_t"]==2 & 
                                   processed_data_sub$Zdata[,"educlevel_t"]==2))] = 2
homophily_edu_f[as.numeric(which(processed_data_sub$Xdata[,"educlevel_t"]==3 & 
                                   processed_data_sub$Zdata[,"educlevel_t"]==3))] = 3
homophily_edu_f[as.numeric(which(processed_data_sub$Xdata[,"educlevel_t"]==4 & 
                                   processed_data_sub$Zdata[,"educlevel_t"]==4))] = 4

# nodematch race
homophily_race_f = matrix(0,nrow=nrow(processed_data_sub$Xdata),ncol=1)
homophily_race_f[as.numeric(which(processed_data_sub$Xdata[,"race_t"]==1 & 
                                   processed_data_sub$Zdata[,"race_t"]==1))] = 1
homophily_race_f[as.numeric(which(processed_data_sub$Xdata[,"race_t"]==2 & 
                                   processed_data_sub$Zdata[,"race_t"]==2))] = 2
homophily_race_f[as.numeric(which(processed_data_sub$Xdata[,"race_t"]==3 & 
                                   processed_data_sub$Zdata[,"race_t"]==3))] = 3
homophily_race_f[as.numeric(which(processed_data_sub$Xdata[,"race_t"]==4 & 
                                   processed_data_sub$Zdata[,"race_t"]==4))] = 4

# nodematch age
homophily_age_f = matrix(0,nrow=nrow(processed_data_sub$Xdata),ncol=1)
homophily_age_f[as.numeric(which(processed_data_sub$Xdata[,"tage_t"]==1 & 
                                    processed_data_sub$Zdata[,"tage_t"]==1))] = 1
homophily_age_f[as.numeric(which(processed_data_sub$Xdata[,"tage_t"]==2 & 
                                    processed_data_sub$Zdata[,"tage_t"]==2))] = 2
homophily_age_f[as.numeric(which(processed_data_sub$Xdata[,"tage_t"]==3 & 
                                    processed_data_sub$Zdata[,"tage_t"]==3))] = 3
homophily_age_f[as.numeric(which(processed_data_sub$Xdata[,"tage_t"]==4 & 
                                    processed_data_sub$Zdata[,"tage_t"]==4))] = 4


females = data.frame(time = diag(processed_data_sub$mu),
                     dissolved = as.numeric(diag(processed_data_sub$mu) < 4),
                     # homophily_less_than_HS = as.numeric(processed_data_sub$Xdata[,"educlevel_t"]==1 & 
                     #                                     processed_data_sub$Zdata[,"educlevel_t"]==1),
                     # homophily_HS = as.numeric(processed_data_sub$Xdata[,"educlevel_t"]==2 & 
                     #                           processed_data_sub$Zdata[,"educlevel_t"]==2),
                     # homophily_some_college = as.numeric(processed_data_sub$Xdata[,"educlevel_t"]==3 & 
                     #                                     processed_data_sub$Zdata[,"educlevel_t"]==3),
                     # homophily_BA_and_beyond = as.numeric(processed_data_sub$Xdata[,"educlevel_t"]==4 & 
                     #                                      processed_data_sub$Zdata[,"educlevel_t"]==4),
                     # 
                     # homophily_Hispanic = as.numeric(processed_data_sub$Xdata[,"race_t"]==1 & 
                     #                                 processed_data_sub$Zdata[,"race_t"]==1),
                     # homophily_Black = as.numeric(processed_data_sub$Xdata[,"race_t"]==2 & 
                     #                              processed_data_sub$Zdata[,"race_t"]==2),
                     # homophily_White = as.numeric(processed_data_sub$Xdata[,"race_t"]==3 & 
                     #                              processed_data_sub$Zdata[,"race_t"]==3),
                     # homophily_Asian = as.numeric(processed_data_sub$Xdata[,"race_t"]==4 & 
                     #                              processed_data_sub$Zdata[,"race_t"]==4),
                     # 
                     # homophily_18_27.9 = as.numeric(processed_data_sub$Xdata[,"tage_t"]==1 & 
                     #                                processed_data_sub$Zdata[,"tage_t"]==1),
                     # homophily_28_38.9 = as.numeric(processed_data_sub$Xdata[,"tage_t"]==2 & 
                     #                                processed_data_sub$Zdata[,"tage_t"]==2),
                     # homophily_39_48.9 = as.numeric(processed_data_sub$Xdata[,"tage_t"]==3 & 
                     #                                processed_data_sub$Zdata[,"tage_t"]==3),
                     # homophily_49_59 = as.numeric(processed_data_sub$Xdata[,"tage_t"]==4 & 
                     #                              processed_data_sub$Zdata[,"tage_t"]==4),
                     # 
                     # greater_than_own_age = as.numeric((processed_data_sub$Xdata[,"tage_t"] ==1 & 
                     #                                    processed_data_sub$Zdata[,"tage_t"] > 1) |
                     #                                   (processed_data_sub$Xdata[,"tage_t"] ==2 & 
                     #                                    processed_data_sub$Zdata[,"tage_t"] > 2) | 
                     #                                   (processed_data_sub$Xdata[,"tage_t"] ==3 & 
                     #                                    processed_data_sub$Zdata[,"tage_t"] > 3)),
                     # homophily_age = as.numeric((processed_data_sub$Xdata[,"tage_t"] == 1 & 
                     #                             processed_data_sub$Zdata[,"tage_t"] == 1) |
                     #                            (processed_data_sub$Xdata[,"tage_t"] == 2 & 
                     #                             processed_data_sub$Zdata[,"tage_t"] == 2) | 
                     #                            (processed_data_sub$Xdata[,"tage_t"] == 3 & 
                     #                             processed_data_sub$Zdata[,"tage_t"] == 3) |
                     #                            (processed_data_sub$Xdata[,"tage_t"] == 4 & 
                     #                             processed_data_sub$Zdata[,"tage_t"] == 4)),
                     # less_than_own_age = as.numeric((processed_data_sub$Xdata[,"tage_t"] ==4 & 
                     #                                 processed_data_sub$Zdata[,"tage_t"] < 4) |
                     #                                (processed_data_sub$Xdata[,"tage_t"] ==3 & 
                     #                                 processed_data_sub$Zdata[,"tage_t"] < 3) | 
                     #                                (processed_data_sub$Xdata[,"tage_t"] ==2 & 
                     #                                 processed_data_sub$Zdata[,"tage_t"] < 2)),
                     
                     age_disparity = age_disparity_f,
                     homophily_edu = homophily_edu_f,
                     homophily_race = homophily_race_f,
                     homophily_age = homophily_age_f
                     )


# sanity check
# cbind(processed_data_sub$Xdata[,"tage_t"], processed_data_sub$Zdata[,"tage_t"])[which(females$less_than_own_age==1),]


# use one variable to indicate age homphily, greater than own, or less than own
age_disparity_m = matrix(0,nrow=nrow(processed_data_sub$Zdata),ncol=1)
idx_greater = as.numeric(which((processed_data_sub$Zdata[,"tage_t"] ==1 & 
                                       processed_data_sub$Xdata[,"tage_t"] > 1) |
                                      (processed_data_sub$Zdata[,"tage_t"] ==2 & 
                                         processed_data_sub$Xdata[,"tage_t"] > 2) | 
                                      (processed_data_sub$Zdata[,"tage_t"] ==3 & 
                                         processed_data_sub$Xdata[,"tage_t"] > 3)))
age_disparity_m[idx_greater] = 1

idx_less = as.numeric(which((processed_data_sub$Zdata[,"tage_t"] ==4 & 
                               processed_data_sub$Xdata[,"tage_t"] < 4) |
                              (processed_data_sub$Zdata[,"tage_t"] ==3 & 
                                 processed_data_sub$Xdata[,"tage_t"] < 3) | 
                              (processed_data_sub$Zdata[,"tage_t"] ==2 & 
                                 processed_data_sub$Xdata[,"tage_t"] < 2)))
age_disparity_m[idx_less] = 2

# nodematch edu
homophily_edu_m = matrix(0,nrow=nrow(processed_data_sub$Zdata),ncol=1)
homophily_edu_m[as.numeric(which(processed_data_sub$Zdata[,"educlevel_t"]==1 & 
                                   processed_data_sub$Xdata[,"educlevel_t"]==1))] = 1
homophily_edu_m[as.numeric(which(processed_data_sub$Zdata[,"educlevel_t"]==2 & 
                                   processed_data_sub$Xdata[,"educlevel_t"]==2))] = 2
homophily_edu_m[as.numeric(which(processed_data_sub$Zdata[,"educlevel_t"]==3 & 
                                   processed_data_sub$Xdata[,"educlevel_t"]==3))] = 3
homophily_edu_m[as.numeric(which(processed_data_sub$Zdata[,"educlevel_t"]==4 & 
                                   processed_data_sub$Xdata[,"educlevel_t"]==4))] = 4

# nodematch race
homophily_race_m = matrix(0,nrow=nrow(processed_data_sub$Zdata),ncol=1)
homophily_race_m[as.numeric(which(processed_data_sub$Zdata[,"race_t"]==1 & 
                                    processed_data_sub$Xdata[,"race_t"]==1))] = 1
homophily_race_m[as.numeric(which(processed_data_sub$Zdata[,"race_t"]==2 & 
                                    processed_data_sub$Xdata[,"race_t"]==2))] = 2
homophily_race_m[as.numeric(which(processed_data_sub$Zdata[,"race_t"]==3 & 
                                    processed_data_sub$Xdata[,"race_t"]==3))] = 3
homophily_race_m[as.numeric(which(processed_data_sub$Zdata[,"race_t"]==4 & 
                                    processed_data_sub$Xdata[,"race_t"]==4))] = 4

# nodematch age
homophily_age_m = matrix(0,nrow=nrow(processed_data_sub$Zdata),ncol=1)
homophily_age_m[as.numeric(which(processed_data_sub$Zdata[,"tage_t"]==1 & 
                                   processed_data_sub$Xdata[,"tage_t"]==1))] = 1
homophily_age_m[as.numeric(which(processed_data_sub$Zdata[,"tage_t"]==2 & 
                                   processed_data_sub$Xdata[,"tage_t"]==2))] = 2
homophily_age_m[as.numeric(which(processed_data_sub$Zdata[,"tage_t"]==3 & 
                                   processed_data_sub$Xdata[,"tage_t"]==3))] = 3
homophily_age_m[as.numeric(which(processed_data_sub$Zdata[,"tage_t"]==4 & 
                                   processed_data_sub$Xdata[,"tage_t"]==4))] = 4

males = data.frame(time = diag(processed_data_sub$mu),
                     dissolved = as.numeric(diag(processed_data_sub$mu) < 4),
                     # homophily_less_than_HS = as.numeric(processed_data_sub$Zdata[,"educlevel_t"]==1 & 
                     #                                       processed_data_sub$Xdata[,"educlevel_t"]==1),
                     # homophily_HS = as.numeric(processed_data_sub$Zdata[,"educlevel_t"]==2 & 
                     #                             processed_data_sub$Xdata[,"educlevel_t"]==2),
                     # homophily_some_college = as.numeric(processed_data_sub$Zdata[,"educlevel_t"]==3 & 
                     #                                       processed_data_sub$Xdata[,"educlevel_t"]==3),
                     # homophily_BA_and_beyond = as.numeric(processed_data_sub$Zdata[,"educlevel_t"]==4 & 
                     #                                        processed_data_sub$Xdata[,"educlevel_t"]==4),
                     # 
                     # homophily_Hispanic = as.numeric(processed_data_sub$Zdata[,"race_t"]==1 & 
                     #                                   processed_data_sub$Xdata[,"race_t"]==1),
                     # homophily_Black = as.numeric(processed_data_sub$Zdata[,"race_t"]==2 & 
                     #                                processed_data_sub$Xdata[,"race_t"]==2),
                     # homophily_White = as.numeric(processed_data_sub$Zdata[,"race_t"]==3 & 
                     #                                processed_data_sub$Xdata[,"race_t"]==3),
                     # homophily_Asian = as.numeric(processed_data_sub$Zdata[,"race_t"]==4 & 
                     #                                processed_data_sub$Xdata[,"race_t"]==4),
                     # 
                     # homophily_18_27.9 = as.numeric(processed_data_sub$Zdata[,"tage_t"]==1 & 
                     #                                  processed_data_sub$Xdata[,"tage_t"]==1),
                     # homophily_28_38.9 = as.numeric(processed_data_sub$Zdata[,"tage_t"]==2 & 
                     #                                  processed_data_sub$Xdata[,"tage_t"]==2),
                     # homophily_39_48.9 = as.numeric(processed_data_sub$Zdata[,"tage_t"]==3 & 
                     #                                  processed_data_sub$Xdata[,"tage_t"]==3),
                     # homophily_49_59 = as.numeric(processed_data_sub$Zdata[,"tage_t"]==4 & 
                     #                                processed_data_sub$Xdata[,"tage_t"]==4),
                     # 
                     # greater_than_own_age = as.numeric((processed_data_sub$Zdata[,"tage_t"] ==1 & 
                     #                                      processed_data_sub$Xdata[,"tage_t"] > 1) |
                     #                                     (processed_data_sub$Zdata[,"tage_t"] ==2 & 
                     #                                        processed_data_sub$Xdata[,"tage_t"] > 2) | 
                     #                                     (processed_data_sub$Zdata[,"tage_t"] ==3 & 
                     #                                        processed_data_sub$Xdata[,"tage_t"] > 3)),
                     # homophily_age = as.numeric((processed_data_sub$Zdata[,"tage_t"] == 1 & 
                     #                               processed_data_sub$Xdata[,"tage_t"] == 1) |
                     #                              (processed_data_sub$Zdata[,"tage_t"] == 2 & 
                     #                                 processed_data_sub$Xdata[,"tage_t"] == 2) | 
                     #                              (processed_data_sub$Zdata[,"tage_t"] == 3 & 
                     #                                 processed_data_sub$Xdata[,"tage_t"] == 3) |
                     #                              (processed_data_sub$Zdata[,"tage_t"] == 4 & 
                     #                                 processed_data_sub$Xdata[,"tage_t"] == 4)),
                     # less_than_own_age = as.numeric((processed_data_sub$Zdata[,"tage_t"] ==4 & 
                     #                                   processed_data_sub$Xdata[,"tage_t"] < 4) |
                     #                                  (processed_data_sub$Zdata[,"tage_t"] ==3 & 
                     #                                     processed_data_sub$Xdata[,"tage_t"] < 3) | 
                     #                                  (processed_data_sub$Zdata[,"tage_t"] ==2 & 
                     #                                     processed_data_sub$Xdata[,"tage_t"] < 2)),
                   
                   age_disparity = age_disparity_m,
                   homophily_edu = homophily_edu_m,
                   homophily_race = homophily_race_m,
                   homophily_age = homophily_age_m
                  )

# females$homophily_less_than_HS = as.factor(females$homophily_less_than_HS)
# females$homophily_HS = as.factor(females$homophily_HS)
# females$homophily_some_college = as.factor(females$homophily_some_college)
# females$homophily_BA_and_beyond = as.factor(females$homophily_BA_and_beyond)
# females$homophily_Hispanic = as.factor(females$homophily_Hispanic)
# females$homophily_Black = as.factor(females$homophily_Black)
# females$homophily_White = as.factor(females$homophily_White)
# females$homophily_Asian = as.factor(females$homophily_Asian)
# females$homophily_18_27.9 = as.factor(females$homophily_18_27.9)
# females$homophily_28_38.9 = as.factor(females$homophily_28_38.9)
# females$homophily_39_48.9 = as.factor(females$homophily_39_48.9)
# females$homophily_49_59 = as.factor(females$homophily_49_59)
# females$greater_than_own_age = as.factor(females$greater_than_own_age)
# females$homophily_age = as.factor(females$homophily_age)
# females$less_than_own_age = as.factor(females$less_than_own_age)
females$age_disparity = as.factor(females$age_disparity)   
females$homophily_edu = as.factor(females$homophily_edu)
females$homophily_race = as.factor(females$homophily_race)
females$homophily_age = as.factor(females$homophily_age)
     
# males$homophily_less_than_HS = as.factor(males$homophily_less_than_HS)
# males$homophily_HS = as.factor(males$homophily_HS)
# males$homophily_some_college = as.factor(males$homophily_some_college)
# males$homophily_BA_and_beyond = as.factor(males$homophily_BA_and_beyond)
# males$homophily_Hispanic = as.factor(males$homophily_Hispanic)
# males$homophily_Black = as.factor(males$homophily_Black)
# males$homophily_White = as.factor(males$homophily_White)
# males$homophily_Asian = as.factor(males$homophily_Asian)
# males$homophily_18_27.9 = as.factor(males$homophily_18_27.9)
# males$homophily_28_38.9 = as.factor(males$homophily_28_38.9)
# males$homophily_39_48.9 = as.factor(males$homophily_39_48.9)
# males$homophily_49_59 = as.factor(males$homophily_49_59)
# males$greater_than_own_age = as.factor(males$greater_than_own_age)
# males$homophily_age = as.factor(males$homophily_age)
# males$less_than_own_age = as.factor(males$less_than_own_age)
males$age_disparity = as.factor(males$age_disparity)          
males$homophily_edu = as.factor(males$homophily_edu)
males$homophily_race = as.factor(males$homophily_race)
males$homophily_age = as.factor(males$homophily_age)

###############################################
#       Cox proportional hazards model        #
###############################################

######### nodematch education #############

edu_surv_object_X <- Surv(time = females$time, event = females$dissolved)
# edu_fit_X <- coxph(edu_surv_object_X ~ homophily_less_than_HS + homophily_HS + homophily_some_college + homophily_BA_and_beyond, 
#                    data = females)
edu_fit_X <- coxph(edu_surv_object_X ~ homophily_edu, data = females)
print(ggforest(edu_fit_X, data = females))
print(summary(edu_fit_X))

# same as fit_X
# edu_surv_object_Y <- Surv(time = males$time, event = males$dissolved)
# edu_fit_Y <- coxph(edu_surv_object_Y ~ homophily_edu, 
#                data = males)
# ggforest(edu_fit_Y, data = males)
# summary(edu_fit_Y)


######### nodematch race #############

race_surv_object_X <- Surv(time = females$time, event = females$dissolved)
# race_fit_X <- coxph(race_surv_object_X ~ homophily_Hispanic + homophily_Black + homophily_White + homophily_Asian, 
#                data = females)
race_fit_X <- coxph(race_surv_object_X ~ homophily_race, data = females)
print(ggforest(race_fit_X, data = females))
print(summary(race_fit_X))


######### nodematch age #############

age_surv_object_X <- Surv(time = females$time, event = females$dissolved)
# age_fit_X <- coxph(age_surv_object_X ~ homophily_18_27.9 + homophily_28_38.9 + homophily_39_48.9 + homophily_49_59, 
#                     data = females)
age_fit_X <- coxph(age_surv_object_X ~ homophily_age, data = females)
print(ggforest(age_fit_X, data = females))
print(summary(age_fit_X))


######### age disparity #############

# can't do this because of "perfect classification" 
# xtabs(~greater_than_own_age + homophily_age + less_than_own_age, data=females)

age_disp_surv_object_X <- Surv(time = females$time, event = females$dissolved)
# age_disp_fit_X <- coxph(age_disp_surv_object_X ~ greater_than_own_age + homophily_age + less_than_own_age,
#                    data = females)
age_disp_fit_X <- coxph(age_disp_surv_object_X ~ age_disparity,
                        data = females)
print(ggforest(age_disp_fit_X, data = females))
print(summary(age_disp_fit_X))


age_disp_surv_object_Y <- Surv(time = males$time, event = males$dissolved)
age_disp_fit_Y <- coxph(age_disp_surv_object_Y ~ age_disparity, data = males)
print(ggforest(age_disp_fit_Y, data = males))
print(summary(age_disp_fit_Y))



###############################################
#          the Kaplan-Meier method            #
###############################################

# I think this method doesn't give reasonable result because everyone
# survived to time = 4 gets censored out at time = 4 (this is my guess)

# # age disparity
# surv_object_X <- Surv(time = females$time, event = females$dissolved) 
# surv_object_Y <- Surv(time = males$time, event = males$dissolved)
# 
# KM_age_disp_X <- survfit(surv_object_X ~ age_disparity, data = females)
# summary(KM_age_disp_X)
# ggsurvplot(KM_age_disp_X, data = females, pval = TRUE)
# 
# # nodematch race
# KM_race_fit_X <- survfit(surv_object_X ~ homophily_Hispanic + homophily_Black + homophily_White + homophily_Asian, 
#                     data = females)
# ggsurvplot(KM_race_fit_X, data = females)
# summary(KM_race_fit_X)

##############################################

# randomly sample the data to show that the estimates would change with availability
# https://stats.stackexchange.com/questions/127471/how-to-make-prediction-in-survival-analysis-using-r
# https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/

source("fitrpm_R_CP.R")
source("loglikelihood_CP.R")
source("equality_constraint_CP.R")
source("rpm.model.matrix.R")
source("choice_probability.R")
source("check_CP.R")
source("asymptotic_var.R")

set.seed(1234)

bootstrap_iter = 20
sample_sizes = c(200, 500, 1000)

# parameters for two-sided revealed preference model
control = list("algorithm"="NLOPT_LD_SLSQP", "symmetric"=FALSE, "sampling_protocol"="INDIV",
               "xtol_rel"=1.0e-8, "print_level"=0,"maxeval"=1000, 
               "ftol_rel"=1.0e-8,"check_derivatives"=FALSE,"ftol_abs"=1.0e-6,"hessian"=FALSE) 
theta_0 = NULL
ff = ~ b1nodematch("race_t") + b2nodematch("race_t")

# race
bootstrap_race_survival = matrix(0, nrow=bootstrap_iter, ncol=4)
bootstrap_race_rpm = matrix(0, nrow=bootstrap_iter, ncol=4)
for (iter in 1:bootstrap_iter) {
  
  sample_size = sample(sample_sizes, 5, replace=T)
  
  idx_0 = which(females$homophily_race==0)
  idx_0 = sample(idx_0, sample_size[1], replace=T)
  idx_1 = which(females$homophily_race==1)
  idx_1 = sample(idx_1, sample_size[2], replace=T)
  idx_2 = which(females$homophily_race==2)
  idx_2 = sample(idx_2, sample_size[3], replace=T)
  idx_3 = which(females$homophily_race==3)
  idx_3 = sample(idx_3, sample_size[4], replace=T)
  idx_4 = which(females$homophily_race==4)
  idx_4 = sample(idx_4, sample_size[5], replace=T)
  
  keep_idx = c(idx_0, idx_1, idx_2, idx_3, idx_4)
  
  # Cox model
  females_sub = females[keep_idx,]
  race_surv_object_X_sub <- Surv(time = females_sub$time, event = females_sub$dissolved)
  race_fit_X_sub <- coxph(race_surv_object_X_sub ~ homophily_race, data = females_sub)
  bootstrap_race_survival[iter,] = race_fit_X_sub$coefficients
  
  # two-sided revealed preference model
  Xdata = processed_data_sub$Xdata[keep_idx,]
  Zdata = processed_data_sub$Zdata[keep_idx,]
  X_w = processed_data_sub$X_w[keep_idx]
  Z_w = processed_data_sub$Z_w[keep_idx]
  # pair_w = processed_datasub$pair_w[keep_idx]
  pair_w = rep(1, length(keep_idx))
  mu = Diagonal(length(keep_idx)) * diag(processed_data_sub$mu)[keep_idx]
  out <- fitrpm_R_CP(ff, mu, Xdata, Zdata, X_w, Z_w, pair_w, theta_0=NULL, choices=4, control=control)
  bootstrap_race_rpm[iter,] = out$solution[c(2,3,4,5)]

}

boostrap_mean_survival = colMeans(bootstrap_race_survival)
boostrap_sd_survival = apply(bootstrap_race_survival, 2, sd)
boostrap_mean_rpm = colMeans(bootstrap_race_rpm)
boostrap_sd_rpm = apply(bootstrap_race_rpm, 2, sd)

# print(ggforest(race_fit_X_sub, data = females_sub))
# print(summary(race_fit_X_sub))



######### histograms of duration of unions in data ############

# education homophily
edu_homophily_table_f = table(females$homophily_edu, females$time)
edu_homophily_table_f = apply(edu_homophily_table_f, 2, function(x) x/rowSums(edu_homophily_table_f))
edu_homophily_table_f_df = data.frame(type = factor(c("reference", "<HS", "HS", "someCollege", "BA+"), 
                                                    levels = c("reference", "<HS", "HS", "someCollege", "BA+")), 
                                      proportion = edu_homophily_table_f[,4])
print(ggplot(edu_homophily_table_f_df, aes(type, proportion)) + geom_bar(stat="identity", position="dodge"))


# race homophily
race_homophily_table_f = table(females$homophily_race, females$time)
race_homophily_table_f = apply(race_homophily_table_f, 2, function(x) x/rowSums(race_homophily_table_f))
race_homophily_table_f_df = data.frame(type = factor(c("reference", "Hispanic", "Black", "White", "Asian"), 
                                                    levels = c("reference", "Hispanic", "Black", "White", "Asian")), 
                                      proportion = race_homophily_table_f[,4])
print(ggplot(race_homophily_table_f_df, aes(type, proportion)) + geom_bar(stat="identity", position="dodge"))


# age homophily
age_homophily_table_f = table(females$homophily_age, females$time)
age_homophily_table_f = apply(age_homophily_table_f, 2, function(x) x/rowSums(age_homophily_table_f))
age_homophily_table_f_df = data.frame(type = factor(c("reference", "18-27.9", "28-38.9", "39-48.9", "49-59"), 
                                                    levels = c("reference", "18-27.9", "28-38.9", "39-48.9", "49-59")), 
                                      proportion = age_homophily_table_f[,4])
print(ggplot(age_homophily_table_f_df, aes(type, proportion)) + geom_bar(stat="identity", position="dodge"))


# female age disparity
age_disp_table_f = table(females$age_disparity, females$time)
age_disp_table_f = apply(age_disp_table_f, 2, function(x) x/rowSums(age_disp_table_f))
age_disp_table_f_df = data.frame(type = factor(c("age homophily", "greater than own", "less than own")), 
                                 proportion = age_disp_table_f[,4])
print(ggplot(age_disp_table_f_df, aes(type, proportion)) + geom_bar(stat="identity", position="dodge"))


# male age disparity
age_disp_table_m = table(males$age_disparity, males$time)
age_disp_table_m = apply(age_disp_table_m, 2, function(x) x/rowSums(age_disp_table_m))
age_disp_table_m_df = data.frame(type = factor(c("age homophily", "greater than own", "less than own")), 
                                 proportion = age_disp_table_m[,4])
print(ggplot(age_disp_table_m_df, aes(type, proportion)) + geom_bar(stat="identity", position="dodge"))

# ggplot(females, aes(time, fill=age_disparity)) + geom_histogram(position="dodge", stat="count")
