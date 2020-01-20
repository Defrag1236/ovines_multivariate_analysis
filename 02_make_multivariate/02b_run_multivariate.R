### run multivariate ###

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("MV_sst_file.Rdata")

# run mv and save results

library(MultiABEL)

multi_weight <- MultiSummary(sst, index=c(1:3, 5:10), type = "outbred")
save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_weight.Rdata", list=c("multi_weight"))
rm(multi_weight)

multi_fat <- MultiSummary(sst, index=c(4, 13:16), type = "outbred")
save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_fat.Rdata", list=c("multi_fat"))
rm(multi_fat)

multi_muscle <- MultiSummary(sst, index=c(17:19), type = "outbred")
save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_muscle.Rdata", list=c("multi_muscle"))
rm(multi_muscle)

multi_concentration_of_some_substancies <- MultiSummary(sst, index=c(22:24, 42:43), type = "outbred")
save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_concentration_of_some_substancies.Rdata", list=c("multi_concentration_of_some_substancies"))
rm(multi_concentration_of_some_substancies)

multi_retail_colour <- MultiSummary(sst, index=c(25:36), type = "outbred")
save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_retail_colour.Rdata", list=c("multi_retail_colour"))
rm(multi_retail_colour)

multi_fresh_colour <- MultiSummary(sst, index=c(37:39), type = "outbred")
save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_fresh_colour.Rdata", list=c("multi_fresh_colour"))
rm(multi_fresh_colour)

multi_concentration_of_some_acids <- MultiSummary(sst, index=c(44:56), type = "outbred")
save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_concentration_of_some_acids.Rdata", list=c("multi_concentration_of_some_acids"))
rm(multi_concentration_of_some_acids)

multi_ph <- MultiSummary(sst, index=c(40:41), type = "outbred")
save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_ph.Rdata", list=c("multi_ph"))
rm(multi_ph)