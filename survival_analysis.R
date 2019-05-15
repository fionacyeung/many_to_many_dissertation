# setwd("C:\\UCLA\\research_projects\\many_to_many_dissertation")
# source("survival_analysis.R")

# see tutorial at https://www.datacamp.com/community/tutorials/survival-analysis-R

library(survival)
library(survminer)

data(ovarian)

ovarian$rx <- factor(ovarian$rx, 
                     levels = c("1", "2"), 
                     labels = c("A", "B"))
ovarian$resid.ds <- factor(ovarian$resid.ds, 
                           levels = c("1", "2"), 
                           labels = c("no", "yes"))
ovarian$ecog.ps <- factor(ovarian$ecog.ps, 
                          levels = c("1", "2"), 
                          labels = c("good", "bad"))

ovarian$age_group = "young"
ovarian$age_group[ovarian$age >= 50] = "old"
ovarian$age_group <- factor(ovarian$age_group)

# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = ovarian$futime, event = ovarian$fustat) # futime = time of death or censored, fustat = censored or not  
# surv_object

fit1 <- survfit(surv_object ~ rx, data = ovarian)
summary(fit1)
ggsurvplot(fit1, data = ovarian, pval = TRUE)

fit2 <- survfit(surv_object ~ resid.ds, data = ovarian)
ggsurvplot(fit2, data = ovarian, pval = TRUE)


# Fit a Cox proportional hazards model
fit.coxph <- coxph(surv_object ~ rx + resid.ds + age_group + ecog.ps, 
                   data = ovarian)
ggforest(fit.coxph, data = ovarian)


########################################################

set.seed(1234)

load("processed_data_sub.Rdata")


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

females = data.frame(time = diag(processed_data_sub$mu),
                     censored = as.numeric(diag(processed_data_sub$mu) >= 4),
                     homophily_less_than_HS = as.numeric(processed_data_sub$Xdata[,"educlevel_t"]==1 & 
                                                         processed_data_sub$Zdata[,"educlevel_t"]==1),
                     homophily_HS = as.numeric(processed_data_sub$Xdata[,"educlevel_t"]==2 & 
                                               processed_data_sub$Zdata[,"educlevel_t"]==2),
                     homophily_some_college = as.numeric(processed_data_sub$Xdata[,"educlevel_t"]==3 & 
                                                         processed_data_sub$Zdata[,"educlevel_t"]==3),
                     homophily_BA_and_beyond = as.numeric(processed_data_sub$Xdata[,"educlevel_t"]==4 & 
                                                          processed_data_sub$Zdata[,"educlevel_t"]==4),
                     
                     homophily_Hispanic = as.numeric(processed_data_sub$Xdata[,"race_t"]==1 & 
                                                     processed_data_sub$Zdata[,"race_t"]==1),
                     homophily_Black = as.numeric(processed_data_sub$Xdata[,"race_t"]==2 & 
                                                  processed_data_sub$Zdata[,"race_t"]==2),
                     homophily_White = as.numeric(processed_data_sub$Xdata[,"race_t"]==3 & 
                                                  processed_data_sub$Zdata[,"race_t"]==3),
                     homophily_Asian = as.numeric(processed_data_sub$Xdata[,"race_t"]==4 & 
                                                  processed_data_sub$Zdata[,"race_t"]==4),
                     
                     homophily_18_27.9 = as.numeric(processed_data_sub$Xdata[,"tage_t"]==1 & 
                                                    processed_data_sub$Zdata[,"tage_t"]==1),
                     homophily_28_38.9 = as.numeric(processed_data_sub$Xdata[,"tage_t"]==2 & 
                                                    processed_data_sub$Zdata[,"tage_t"]==2),
                     homophily_39_48.9 = as.numeric(processed_data_sub$Xdata[,"tage_t"]==3 & 
                                                    processed_data_sub$Zdata[,"tage_t"]==3),
                     homophily_49_59 = as.numeric(processed_data_sub$Xdata[,"tage_t"]==4 & 
                                                  processed_data_sub$Zdata[,"tage_t"]==4),
                     
                     greater_than_own_age = as.numeric((processed_data_sub$Xdata[,"tage_t"] ==1 & 
                                                        processed_data_sub$Zdata[,"tage_t"] > 1) |
                                                       (processed_data_sub$Xdata[,"tage_t"] ==2 & 
                                                        processed_data_sub$Zdata[,"tage_t"] > 2) | 
                                                       (processed_data_sub$Xdata[,"tage_t"] ==3 & 
                                                        processed_data_sub$Zdata[,"tage_t"] > 3)),
                     homophily_age = as.numeric((processed_data_sub$Xdata[,"tage_t"] == 1 & 
                                                 processed_data_sub$Zdata[,"tage_t"] == 1) |
                                                (processed_data_sub$Xdata[,"tage_t"] == 2 & 
                                                 processed_data_sub$Zdata[,"tage_t"] == 2) | 
                                                (processed_data_sub$Xdata[,"tage_t"] == 3 & 
                                                 processed_data_sub$Zdata[,"tage_t"] == 3) |
                                                (processed_data_sub$Xdata[,"tage_t"] == 4 & 
                                                 processed_data_sub$Zdata[,"tage_t"] == 4)),
                     less_than_own_age = as.numeric((processed_data_sub$Xdata[,"tage_t"] ==4 & 
                                                     processed_data_sub$Zdata[,"tage_t"] < 4) |
                                                    (processed_data_sub$Xdata[,"tage_t"] ==3 & 
                                                     processed_data_sub$Zdata[,"tage_t"] < 3) | 
                                                    (processed_data_sub$Xdata[,"tage_t"] ==2 & 
                                                     processed_data_sub$Zdata[,"tage_t"] < 2)),
                     
                     age_disparity = age_disparity_f
                     )


# sanity check
# cbind(processed_data_sub$Xdata[,"tage_t"], processed_data_sub$Zdata[,"tage_t"])[which(females$less_than_own_age==1),]

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

males = data.frame(time = diag(processed_data_sub$mu),
                     censored = as.numeric(diag(processed_data_sub$mu) >= 4),
                     homophily_less_than_HS = as.numeric(processed_data_sub$Zdata[,"educlevel_t"]==1 & 
                                                           processed_data_sub$Xdata[,"educlevel_t"]==1),
                     homophily_HS = as.numeric(processed_data_sub$Zdata[,"educlevel_t"]==2 & 
                                                 processed_data_sub$Xdata[,"educlevel_t"]==2),
                     homophily_some_college = as.numeric(processed_data_sub$Zdata[,"educlevel_t"]==3 & 
                                                           processed_data_sub$Xdata[,"educlevel_t"]==3),
                     homophily_BA_and_beyond = as.numeric(processed_data_sub$Zdata[,"educlevel_t"]==4 & 
                                                            processed_data_sub$Xdata[,"educlevel_t"]==4),
                     
                     homophily_Hispanic = as.numeric(processed_data_sub$Zdata[,"race_t"]==1 & 
                                                       processed_data_sub$Xdata[,"race_t"]==1),
                     homophily_Black = as.numeric(processed_data_sub$Zdata[,"race_t"]==2 & 
                                                    processed_data_sub$Xdata[,"race_t"]==2),
                     homophily_White = as.numeric(processed_data_sub$Zdata[,"race_t"]==3 & 
                                                    processed_data_sub$Xdata[,"race_t"]==3),
                     homophily_Asian = as.numeric(processed_data_sub$Zdata[,"race_t"]==4 & 
                                                    processed_data_sub$Xdata[,"race_t"]==4),
                     
                     homophily_18_27.9 = as.numeric(processed_data_sub$Zdata[,"tage_t"]==1 & 
                                                      processed_data_sub$Xdata[,"tage_t"]==1),
                     homophily_28_38.9 = as.numeric(processed_data_sub$Zdata[,"tage_t"]==2 & 
                                                      processed_data_sub$Xdata[,"tage_t"]==2),
                     homophily_39_48.9 = as.numeric(processed_data_sub$Zdata[,"tage_t"]==3 & 
                                                      processed_data_sub$Xdata[,"tage_t"]==3),
                     homophily_49_59 = as.numeric(processed_data_sub$Zdata[,"tage_t"]==4 & 
                                                    processed_data_sub$Xdata[,"tage_t"]==4),
                     
                     greater_than_own_age = as.numeric((processed_data_sub$Zdata[,"tage_t"] ==1 & 
                                                          processed_data_sub$Xdata[,"tage_t"] > 1) |
                                                         (processed_data_sub$Zdata[,"tage_t"] ==2 & 
                                                            processed_data_sub$Xdata[,"tage_t"] > 2) | 
                                                         (processed_data_sub$Zdata[,"tage_t"] ==3 & 
                                                            processed_data_sub$Xdata[,"tage_t"] > 3)),
                     homophily_age = as.numeric((processed_data_sub$Zdata[,"tage_t"] == 1 & 
                                                   processed_data_sub$Xdata[,"tage_t"] == 1) |
                                                  (processed_data_sub$Zdata[,"tage_t"] == 2 & 
                                                     processed_data_sub$Xdata[,"tage_t"] == 2) | 
                                                  (processed_data_sub$Zdata[,"tage_t"] == 3 & 
                                                     processed_data_sub$Xdata[,"tage_t"] == 3) |
                                                  (processed_data_sub$Zdata[,"tage_t"] == 4 & 
                                                     processed_data_sub$Xdata[,"tage_t"] == 4)),
                     less_than_own_age = as.numeric((processed_data_sub$Zdata[,"tage_t"] ==4 & 
                                                       processed_data_sub$Xdata[,"tage_t"] < 4) |
                                                      (processed_data_sub$Zdata[,"tage_t"] ==3 & 
                                                         processed_data_sub$Xdata[,"tage_t"] < 3) | 
                                                      (processed_data_sub$Zdata[,"tage_t"] ==2 & 
                                                         processed_data_sub$Xdata[,"tage_t"] < 2)),
                   
                   age_disparity = age_disparity_m
                  )

females$homophily_less_than_HS = as.factor(females$homophily_less_than_HS)
females$homophily_HS = as.factor(females$homophily_HS)
females$homophily_some_college = as.factor(females$homophily_some_college)
females$homophily_BA_and_beyond = as.factor(females$homophily_BA_and_beyond)
females$homophily_Hispanic = as.factor(females$homophily_Hispanic)
females$homophily_Black = as.factor(females$homophily_Black)
females$homophily_White = as.factor(females$homophily_White)
females$homophily_Asian = as.factor(females$homophily_Asian)
females$homophily_18_27.9 = as.factor(females$homophily_18_27.9)
females$homophily_28_38.9 = as.factor(females$homophily_28_38.9)
females$homophily_39_48.9 = as.factor(females$homophily_39_48.9)
females$homophily_49_59 = as.factor(females$homophily_49_59)
females$greater_than_own_age = as.factor(females$greater_than_own_age)
females$homophily_age = as.factor(females$homophily_age)
females$less_than_own_age = as.factor(females$less_than_own_age)
females$age_disparity = as.factor(females$age_disparity)             
     
males$homophily_less_than_HS = as.factor(males$homophily_less_than_HS)
males$homophily_HS = as.factor(males$homophily_HS)
males$homophily_some_college = as.factor(males$homophily_some_college)
males$homophily_BA_and_beyond = as.factor(males$homophily_BA_and_beyond)
males$homophily_Hispanic = as.factor(males$homophily_Hispanic)
males$homophily_Black = as.factor(males$homophily_Black)
males$homophily_White = as.factor(males$homophily_White)
males$homophily_Asian = as.factor(males$homophily_Asian)
males$homophily_18_27.9 = as.factor(males$homophily_18_27.9)
males$homophily_28_38.9 = as.factor(males$homophily_28_38.9)
males$homophily_39_48.9 = as.factor(males$homophily_39_48.9)
males$homophily_49_59 = as.factor(males$homophily_49_59)
males$greater_than_own_age = as.factor(males$greater_than_own_age)
males$homophily_age = as.factor(males$homophily_age)
males$less_than_own_age = as.factor(males$less_than_own_age)
males$age_disparity = as.factor(males$age_disparity)          

######### nodematch education #############

edu_surv_object_X <- Surv(time = females$time, event = females$censored)
edu_fit_X <- coxph(edu_surv_object_X ~ homophily_less_than_HS + homophily_HS + homophily_some_college + homophily_BA_and_beyond, 
                   data = females)
ggforest(edu_fit_X, data = females)
summary(edu_fit_X)

# same as fit_X
# edu_surv_object_Y <- Surv(time = males$time, event = males$censored)
# edu_fit_Y <- coxph(edu_surv_object_Y ~ homophily_less_than_HS + homophily_HS + homophily_some_college + homophily_BA_and_beyond, 
#                data = males)
# ggforest(edu_fit_Y, data = males)
# summary(edu_fit_Y)


######### nodematch race #############

race_surv_object_X <- Surv(time = females$time, event = females$censored)
race_fit_X <- coxph(race_surv_object_X ~ homophily_Hispanic + homophily_Black + homophily_White + homophily_Asian, 
               data = females)
ggforest(race_fit_X, data = females)
summary(race_fit_X)


######### nodematch age #############

age_surv_object_X <- Surv(time = females$time, event = females$censored)
age_fit_X <- coxph(age_surv_object_X ~ homophily_18_27.9 + homophily_28_38.9 + homophily_39_48.9 + homophily_49_59, 
                    data = females)
ggforest(age_fit_X, data = females)
summary(age_fit_X)


######### age disparity #############

# can't do this because of "perfect classification" 
# xtabs(~greater_than_own_age + homophily_age + less_than_own_age, data=females)

age_disp_surv_object_X <- Surv(time = females$time, event = females$censored)
# age_disp_fit_X <- coxph(age_disp_surv_object_X ~ greater_than_own_age + homophily_age + less_than_own_age,
#                    data = females)
age_disp_fit_X <- coxph(age_disp_surv_object_X ~ age_disparity,
                        data = females)
ggforest(age_disp_fit_X, data = females)
summary(age_disp_fit_X)


age_disp_surv_object_Y <- Surv(time = males$time, event = males$censored)
age_disp_fit_Y <- coxph(age_disp_surv_object_Y ~ age_disparity, data = males)
ggforest(age_disp_fit_Y, data = males)
summary(age_disp_fit_Y)

