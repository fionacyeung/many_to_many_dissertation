setwd("C:\\UCLA\\research_projects\\many_to_many_dissertation")

library(knitr)
library(ggplot2)

################## nodematch education ####################
rm(list=ls())

load("SE_B1000_nodematch_educlevel_t.RData")

nparam = 5
edu_coef = out$solution[1:(2*nparam)]
edu_coef_f = edu_coef[1:nparam]
edu_coef_m = edu_coef[nparam + 1:nparam]
se_f = se[1:nparam]
se_m = se[nparam+(1:nparam)]

edu_df = data.frame(coef = c(edu_coef_f[-1], edu_coef_m[-1]), name=rep(c("<HS", "HS", "SomeCollege", "BA+"),2),
                    gender=rep(c("female", "male"), each=nparam-1), 
                    se = c(se_f[-1], se_m[-1]))

# output table
print("estimates and bootstrap standard error:")
dff = data.frame(cbind(out$solution, se, sqrt(out$covar2)))
print(kable(dff, col.names = c("estimates", "bootstrap SE", "asympt. SE")))

print(ggplot(edu_df, aes(x=name, y=coef, fill=gender)) + geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(limits=c("<HS", "HS", "SomeCollege", "BA+")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("Homophily") +
  geom_errorbar(aes(ymin=coef-se, ymax=coef+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)))


# plot data distribution
edu_pairs_df_sub = data.frame(pair_id = processed_data_sub$females[,"pair_id"], female_edu = as.factor(processed_data_sub$females[,"educlevel_t"]), 
                              male_edu = as.factor(processed_data_sub$males[,"educlevel_t"]), type = rep(NA, nrow(processed_data_sub$females))) 

# Rename all levels
levels(edu_pairs_df_sub$female_edu) <- paste0("F:",c("<HS", "HS", "SomeCollege", "BA+"))
levels(edu_pairs_df_sub$male_edu) <- paste0("M:",c("<HS", "HS", "SomeCollege", "BA+"))
ggplot(edu_pairs_df_sub, aes(x=male_edu)) + geom_histogram(position="dodge", stat="count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ female_edu) + xlab("")

# heatmap
edu_table = table(edu_pairs_df_sub$female_edu, edu_pairs_df_sub$male_edu)
edu_table_df = data.frame(female_edu=rep(c("<HS", "HS", "SomeCollege", "BA+"), times=4), 
                          male_edu=rep(c("<HS", "HS", "SomeCollege", "BA+"), each=4), 
                          count=matrix(edu_table, ncol = 1))
edu_table_df$female_edu = factor(edu_table_df$female_edu, levels=c("<HS", "HS", "SomeCollege", "BA+"))
edu_table_df$male_edu = factor(edu_table_df$male_edu, levels=c("<HS", "HS", "SomeCollege", "BA+"))

ggplot(edu_table_df, aes(male_edu, female_edu)) + geom_tile(aes(fill = count), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Male education") + ylab("Female education") + geom_text(aes(male_edu, female_edu, label = count), color = "black", size = 4) 


################## nodematch race ########################

rm(list=ls())

load("SE_B1000_nodematch_race_t.RData")

nparam = 5
race_coef = out$solution[1:(2*nparam)]
race_coef_f = race_coef[1:nparam]
race_coef_m = race_coef[nparam + 1:nparam]
se_f = se[1:nparam]
se_m = se[nparam+(1:nparam)]

race_df = data.frame(coef = c(race_coef_f[-1], race_coef_m[-1]), name=rep(c("Hispanic", "Black", "White", "Asian"),2),
                     gender=rep(c("female", "mMale"), each=nparam-1),
                     se = c(se_f[-1], se_m[-1]))

# output table
print("estimates and bootstrap standard error:")
dff = data.frame(cbind(out$solution, se, sqrt(out$covar2)))
print(kable(dff, col.names = c("estimates", "bootstrap SE", "asympt. SE")))

print(ggplot(race_df, aes(x=name, y=coef, fill=gender)) + geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(limits=c("Hispanic", "Black", "White", "Asian")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("Homophily") +
  geom_errorbar(aes(ymin=coef-se, ymax=coef+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)))


# plot the data distribution
race_pairs_df_sub = data.frame(pair_id = processed_data_sub$females[,"pair_id"], female_race = as.factor(processed_data_sub$females[,"race_t"]), 
                               male_race = as.factor(processed_data_sub$males[,"race_t"]), type = rep(NA, nrow(processed_data_sub$females))) 

# Rename all levels
levels(race_pairs_df_sub$female_race) <- paste0("F:",c("Hispanic", "Black", "White", "Asian"))
levels(race_pairs_df_sub$male_race) <- paste0("M:",c("Hispanic", "Black", "White", "Asian"))
ggplot(race_pairs_df_sub, aes(x=male_race)) + geom_histogram(position="dodge", stat="count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ female_race) + xlab("")

# heatmap
race_table = table(race_pairs_df_sub$female_race, race_pairs_df_sub$male_race)
race_table_df = data.frame(female_race=rep(c("Hispanic", "Black", "White", "Asian"), times=4), 
                          male_race=rep(c("Hispanic", "Black", "White", "Asian"), each=4), 
                          count=matrix(race_table, ncol = 1))
race_table_df$female_race = factor(race_table_df$female_race, levels=c("Hispanic", "Black", "White", "Asian"))
race_table_df$male_race = factor(race_table_df$male_race, levels=c("Hispanic", "Black", "White", "Asian"))

ggplot(race_table_df, aes(male_race, female_race)) + geom_tile(aes(fill = count), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Male race") + ylab("Female race") + geom_text(aes(male_race, female_race, label = count), color = "black", size = 4) 


################## nodematch age ########################

rm(list=ls())

load("SE_B1000_nodematch_tage_t.RData")

nparam = 5
age_coef = out$solution[1:(2*nparam)]
age_coef_f = age_coef[1:nparam]
age_coef_m = age_coef[nparam + 1:nparam]
se_f = se[1:nparam]
se_m = se[nparam+(1:nparam)]

age_df = data.frame(coef = c(age_coef_f[-1], age_coef_m[-1]), name=rep(c("18-27.9", "28-38.9", "39-48.9", "49-59"),2),
                    gender=rep(c("female", "male"), each=nparam-1),
                    se = c(se_f[-1], se_m[-1]))

# output table
print("estimates and bootstrap standard error:")
dff = data.frame(cbind(out$solution, se, sqrt(out$covar2)))
print(kable(dff, col.names = c("estimates", "bootstrap SE", "asympt. SE")))

print(ggplot(age_df, aes(x=name, y=coef, fill=gender)) + geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(limits=c("18-27.9", "28-38.9", "39-48.9", "49-59")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("Homophily") +
  geom_errorbar(aes(ymin=coef-se, ymax=coef+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)))


# plot data distribution
age_pairs_df_sub = data.frame(pair_id = processed_data_sub$females[,"pair_id"], female_age = processed_data_sub$females[,"tage_t"], 
                              male_age = processed_data_sub$males[,"tage_t"]) 

age_pairs_df_sub$female_age = as.factor(age_pairs_df_sub$female_age)
age_pairs_df_sub$male_age = as.factor(age_pairs_df_sub$male_age)

# Rename all levels (doesn't work for)
levels(age_pairs_df_sub$female_age) <- c(paste0("F:",c("18-27.9", "28-38.9", "39-48.9", "49-59")))
levels(age_pairs_df_sub$male_age) <- c(paste0("M:",c("18-27.9", "28-38.9", "39-48.9", "49-59")))
ggplot(age_pairs_df_sub, aes(x=male_age)) + geom_histogram(position="dodge", stat="count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ female_age) + xlab("")


# heatmap
age_table = table(age_pairs_df_sub$female_age, age_pairs_df_sub$male_age)
age_table_df = data.frame(female_age=rep(c("18-27.9", "28-38.9", "39-48.9", "49-59"), times=4), 
                           male_age=rep(c("18-27.9", "28-38.9", "39-48.9", "49-59"), each=4), 
                           count=matrix(age_table, ncol = 1))
age_table_df$female_age = factor(age_table_df$female_age, levels=c("18-27.9", "28-38.9", "39-48.9", "49-59"))
age_table_df$male_age = factor(age_table_df$male_age, levels=c("18-27.9", "28-38.9", "39-48.9", "49-59"))

ggplot(age_table_df, aes(male_age, female_age)) + geom_tile(aes(fill = count), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Male age") + ylab("Female age") + geom_text(aes(male_age, female_age, label = count), color = "black", size = 4)


############ age homophily, greater than, smaller than ##############

rm(list=ls())

load("SE_B1000_age_homophily_weight.RData")

nparam = 4
age_coef = out$solution[1:(2*nparam)]
coef_f = age_coef[1:nparam]
coef_m = age_coef[nparam + 1:nparam]
se_f = se[1:nparam]
se_m = se[nparam+(1:nparam)]

age_df = data.frame(coef = c(coef_f[-1], coef_m[-1]),
                    name = rep(c("homophily", 
                                 "greater than own", 
                                 "less than own"), 2),
                    gender = c(rep("female", nparam-1), rep("male", nparam-1)),
                    se = c(se_f[-1], se_m[-1]))

# output table
print("estimates and bootstrap standard error:")
dff = data.frame(cbind(out$solution, se, sqrt(out$covar2)))
print(kable(dff, col.names = c("estimates", "bootstrap SE", "asympt. SE")))

print(ggplot(age_df, aes(x=name, y=coef, fill=gender)) + geom_bar(stat="identity", position=position_dodge(width = 0.9)) + # position="dodge") +
  # scale_x_discrete(limits=age_df$name) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Age") + ylab("Preference coefficient") +
  geom_errorbar(aes(ymin=coef-se, ymax=coef+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)))  


# plot data distribution
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

age_pairs_df_3 = data.frame(pair_type = c(age_pairs_df_2$female_pair_type, age_pairs_df_2$male_pair_type),
                            gender = rep(c("female", "male"), each=nrow(age_pairs_df_2)))

ggplot(age_pairs_df_3, aes(x=pair_type)) + geom_histogram(position="dodge", stat="count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(. ~ gender) + xlab("")

# age_pairs_df_3 = data.frame(pair_type = c(age_pairs_df_2$female_pair_type, age_pairs_df_2$male_pair_type),
#                             gender = rep(c("female", "male"), each=nrow(age_pairs_df_2)),
#                             sample_w = c(processed_data_sub$X_w, processed_data_sub$Z_w))
# 
# ggplot(age_pairs_df_3, aes(x=pair_type)) + geom_histogram(aes(weight = sample_w), position="dodge", stat="count") + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(. ~ gender) + xlab("")
# 

