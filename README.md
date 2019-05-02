# many_to_many_dissertation

This version handles correlated choices. 

unpartnered_partnered_3.R reads in the data and calls fitrpm_R_CP() to do the estimation.

sample_output.txt shows an example of what's printed on the screen when you source unpartnered_partnered_3.R.

processed_data_sub_age_homophily_greater_smaller.Rdata is the final data subset with 1435 pairs.

The variable "B" in unpartnered_partnered_3.R sets the number of iterations used to calculate the standard error of the estimates by bootstrapping.