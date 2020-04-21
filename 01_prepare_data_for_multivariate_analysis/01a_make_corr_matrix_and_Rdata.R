### make correlation matrix and .Rdata file for next steps ###

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/supplement_from_Bolormaa_et_al_2016")

library(data.table)

s4 <- fread("12864_2016_2538_MOESM4_ESM", head=T, stringsAsFactors=F, data.table=F)
s5 <- fread("12864_2016_2538_MOESM5_ESM", head=T, stringsAsFactors=F, data.table=F)
s6 <- fread("12864_2016_2538_MOESM6_ESM", head=T, stringsAsFactors=F, data.table=F)
s7 <- fread("12864_2016_2538_MOESM7_ESM", head=T, stringsAsFactors=F, data.table=F)

# rbind data in 1 data.frame

s_all <- rbind(s4, s5, s6, s7)
rownames(s_all) <- s_all [,1]
s_all <- s_all[,-c(1)]

# make correlation matrix plot for all traits

cor_s_all <- cor(s_all)

library(corrplot)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/plots")

png("correlation_matrix_for_traits_from_Bolormaa.png")

corrplot(cor_s_all)

dev.off()

# make correlation matrix plot for each mv_trait

	## multi_weight

	cor_s_multi_weight <- cor(s_all[,c("PSWT", "HCWT", "LLWT", "TOP", "RND", "BONE", "LEGBONE", "LMY", "DRESSING")])

	png("correlation_matrix_for_multi_weight.png")

	corrplot(cor_s_multi_weight, title="multi_weight", mar=c(0,0,1,0))

	dev.off()

	## multi_fat

	cor_s_multi_fat <- cor(s_all[,c("IMF", "LLFAT", "CFATSCORE", "CCFAT", "HGRFAT", "CFAT5")])

	png("correlation_matrix_for_multi_fat.png")

	corrplot(cor_s_multi_fat, title="multi_fat", mar=c(0,0,1,0))

	dev.off()

	## multi_muscle

	cor_s_multi_muscle <- cor(s_all[,c("CEMW", "CEMD", "CEMA")])

	png("correlation_matrix_for_multi_muscle.png")

	corrplot(cor_s_multi_muscle, title="multi_muscle", mar=c(0,0,1,0))

	dev.off()

	## multi_concentration_of_some_substancies

	cor_s_multi_concentration_of_some_substancies <- cor(s_all[,c("MYOGLOBIN", "GLYCOGEN", "ICDHACTIVITY", "IRONWET", "ZINCWET")])

	png("correlation_matrix_for_multi_concentration_of_some_substancies.png")

	corrplot(cor_s_multi_concentration_of_some_substancies, title="multi_concentration_of_some_substancies", mar=c(0,0,1,0))

	dev.off()

	## multi_retail_colour

	cor_s_multi_retail_colour <- cor(s_all[,c("RCL4", "RCL3", "RCL2", "RCL1", "RCb4", "RCb3", "RCb2", "RCb1", "RCa4", "RCa3", "RCa2", "RCa1")])

	png("correlation_matrix_for_multi_retail_colour.png")

	corrplot(cor_s_multi_retail_colour, title="multi_retail_colour", mar=c(0,0,1,0))

	dev.off()

	## multi_fresh_colour

	cor_s_multi_fresh_colour <- cor(s_all[,c("CFb", "CFa", "CFL")])

	png("correlation_matrix_for_multi_fresh_colour.png")

	corrplot(cor_s_multi_fresh_colour, title="multi_fresh_colour", mar=c(0,0,1,0))

	dev.off()

	## multi_concentration_of_some_acids

	cor_s_multi_concentration_of_some_acids <- cor(s_all[,c("EPADPADHA", "EPADHA", "FA_C22_6n3", "FA_C22_5n3", "FA_C20_5n3", "FA_C20_4n6", "FA_C20_3n6", "FA_C18_2n6", "FA_C18_0", "FA_C16_0", "FA_C14_0", "FA_C12_0", "FA_C10_0")])

	png("correlation_matrix_for_multi_concentration_of_some_acids.png")

	corrplot(cor_s_multi_concentration_of_some_acids, title="multi_concentration_of_some_acids", mar=c(0,0,1,0))

	dev.off()


	## multi_ph

	cor_s_multi_ph <- cor(s_all[,c("PH24LL", "PH24ST")])

	png("correlation_matrix_for_multi_ph.png")

	corrplot(cor_s_multi_ph, title="multi_ph", mar=c(0,0,1,0))

	dev.off()



# save common data.frame as .Rdata

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

save(s_all, file="supplement_from_Bolormaa.Rdata")