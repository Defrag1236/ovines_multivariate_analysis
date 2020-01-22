### rum moltivariate with formela from stephens 2013 ###

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("locus_table_from_Bolormaa.Rdata")

load("supplement_from_Bolormaa.Rdata")

# compare data

s_all_clear <- s_all[rownames(s_all) %in% rownames(locus_table),]
s_all_clear <- s_all_clear[order(rownames(s_all_clear)),]

# run multivariate for all traits

## make corr matrix

cor_s_all <- cor(s_all_clear)

##1st trait (weight)

weight_multivariate <- cbind(locus_table[,c("SNP", "CHR", "POS")], s_all_clear[,c("PSWT", "HCWT", "LLWT", "TOP", "RND", "BONE", "LEGBONE", "LMY", "DRESSING")])

chi_weight_multivariate <- matrix(ncol=1, nrow=nrow(weight_multivariate))
rownames(chi_weight_multivariate) <- rownames(locus_table)

for (n  in (1:nrow(weight_multivariate))) {

	Z <- as.numeric(weight_multivariate[n,-c(1:3)])
	M <- cor_s_all[colnames (weight_multivariate[-c(1:3)]), colnames (weight_multivariate[-c(1:3)])]

	chi_weight_multivariate [n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_weight_multivariate <- matrix(ncol=1, nrow=nrow(weight_multivariate))
p_value_weight_multivariate[,1] <- pchisq(chi_weight_multivariate, df=7, low=F)
rownames(p_value_weight_multivariate) <- rownames(locus_table)

save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_weight.Rdata", list=c("p_value_weight_multivariate"))

##2nd trait (fat)

fat_multivariate <- cbind(locus_table[,c("SNP", "CHR", "POS")], s_all_clear[,c("IMF", "LLFAT", "CFATSCORE", "CCFAT", "HGRFAT", "CFAT5")])

chi_fat_multivariate <- matrix(ncol=1, nrow=nrow(fat_multivariate))
rownames(chi_fat_multivariate) <- rownames(locus_table)

for (n  in (1:nrow(fat_multivariate))) {

	Z <- as.numeric(fat_multivariate[n,-c(1:3)])
	M <- cor_s_all[colnames (fat_multivariate[-c(1:3)]), colnames (fat_multivariate[-c(1:3)])]

	chi_fat_multivariate [n] <- t(Z) %*% solve(M) %*% Z

	}


p_value_fat_multivariate <- matrix(ncol=1, nrow=nrow(fat_multivariate))
p_value_fat_multivariate <- pchisq(chi_fat_multivariate, df=6, low=F)
rownames(p_value_fat_multivariate) <- rownames(locus_table)

save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_fat.Rdata", list=c("p_value_fat_multivariate"))

##3rd trait (muscle)

muscle_multivariate <- cbind(locus_table[,c("SNP", "CHR", "POS")], s_all_clear[,c("CEMW", "CEMD", "CEMA")])

chi_muscle_multivariate <- matrix(ncol=1, nrow=nrow(muscle_multivariate))
rownames(chi_muscle_multivariate) <- rownames(locus_table)

for (n  in (1:nrow(muscle_multivariate))) {

	Z <- as.numeric(muscle_multivariate[n,-c(1:3)])
	M <- cor_s_all[colnames (muscle_multivariate[-c(1:3)]), colnames (muscle_multivariate[-c(1:3)])]

	chi_muscle_multivariate [n] <- t(Z) %*% solve(M) %*% Z

	}


p_value_muscle_multivariate <- matrix(ncol=1, nrow=nrow(muscle_multivariate))
p_value_muscle_multivariate <- pchisq(chi_muscle_multivariate, df=3, low=F)
rownames(p_value_muscle_multivariate) <- rownames(locus_table)

save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_muscle.Rdata", list=c("p_value_muscle_multivariate"))

##4th trait (Concentration of some substancies)

concentration_of_some_substancies_multivariate <- cbind(locus_table[,c("SNP", "CHR", "POS")], s_all_clear[,c("MYOGLOBIN", "GLYCOGEN", "ICDHACTIVITY", "IRONWET", "ZINCWET")])

chi_concentration_of_some_substancies_multivariate <- matrix(ncol=1, nrow=nrow(concentration_of_some_substancies_multivariate))
rownames(chi_concentration_of_some_substancies_multivariate) <- rownames(locus_table)

for (n  in (1:nrow(concentration_of_some_substancies_multivariate))) {

    Z <- as.numeric(concentration_of_some_substancies_multivariate[n,-c(1:3)])
    M <- cor_s_all[colnames (concentration_of_some_substancies_multivariate[-c(1:3)]), colnames (concentration_of_some_substancies_multivariate[-c(1:3)])]

    chi_concentration_of_some_substancies_multivariate [n] <- t(Z) %*% solve(M) %*% Z

    }


p_value_concentration_of_some_substancies_multivariate <- matrix(ncol=1, nrow=nrow(concentration_of_some_substancies_multivariate))
p_value_concentration_of_some_substancies_multivariate <- pchisq(chi_concentration_of_some_substancies_multivariate, df=5, low=F)
rownames(p_value_concentration_of_some_substancies_multivariate) <- rownames(locus_table)

save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_concentration_of_some_substancies.Rdata", list=c("p_value_concentration_of_some_substancies_multivariate"))

##5th trait (retail_colour)

retail_colour_multivariate <- cbind(locus_table[,c("SNP", "CHR", "POS")], s_all_clear[,c("RCL4", "RCL3", "RCL2", "RCL1", "RCb4", "RCb3", "RCb2", "RCb1", "RCa4", "RCa3", "RCa2", "RCa1")])

chi_retail_colour_multivariate <- matrix(ncol=1, nrow=nrow(retail_colour_multivariate))
rownames(chi_retail_colour_multivariate) <- rownames(locus_table)

for (n  in (1:nrow(retail_colour_multivariate))) {

    Z <- as.numeric(retail_colour_multivariate[n,-c(1:3)])
    M <- cor_s_all[colnames (retail_colour_multivariate[-c(1:3)]), colnames (retail_colour_multivariate[-c(1:3)])]

    chi_retail_colour_multivariate [n] <- t(Z) %*% solve(M) %*% Z

    }


p_value_retail_colour_multivariate <- matrix(ncol=1, nrow=nrow(retail_colour_multivariate))
p_value_retail_colour_multivariate <- pchisq(chi_retail_colour_multivariate, df=12, low=F)
rownames(p_value_retail_colour_multivariate) <- rownames(locus_table)

save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_retail_colour.Rdata", list=c("p_value_retail_colour_multivariate"))

##6th trait (fresh_colour)

fresh_colour_multivariate <- cbind(locus_table[,c("SNP", "CHR", "POS")], s_all_clear[,c("CFb", "CFa", "CFL")])

chi_fresh_colour_multivariate <- matrix(ncol=1, nrow=nrow(fresh_colour_multivariate))
rownames(chi_fresh_colour_multivariate) <- rownames(locus_table)

for (n  in (1:nrow(fresh_colour_multivariate))) {

    Z <- as.numeric(fresh_colour_multivariate[n,4:6])
    M <- cor_s_all[colnames (fresh_colour_multivariate[4:6]), colnames (fresh_colour_multivariate[4:6])]

    chi_fresh_colour_multivariate [n] <- t(Z) %*% solve(M) %*% Z

    }


p_value_fresh_colour_multivariate <- matrix(ncol=1, nrow=nrow(fresh_colour_multivariate))
p_value_fresh_colour_multivariate <- pchisq(chi_fresh_colour_multivariate, df=3, low=F)
rownames(p_value_fresh_colour_multivariate) <- rownames(locus_table)

save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_fresh_colour.Rdata", list=c("p_value_fresh_colour_multivariate"))

##7th trait (concentration_of_some_acids)

concentration_of_some_acids_multivariate <- cbind(locus_table[,c("SNP", "CHR", "POS")], s_all_clear[,c("EPADPADHA", "EPADHA", "FA_C22_6n3", "FA_C22_5n3", "FA_C20_5n3", "FA_C20_4n6", "FA_C20_3n6", "FA_C18_2n6", "FA_C18_0", "FA_C16_0", "FA_C14_0", "FA_C12_0", "FA_C10_0")])

chi_concentration_of_some_acids_multivariate <- matrix(ncol=1, nrow=nrow(concentration_of_some_acids_multivariate))
rownames(chi_concentration_of_some_acids_multivariate) <- rownames(locus_table)

for (n  in (1:nrow(concentration_of_some_acids_multivariate))) {

    Z <- as.numeric(concentration_of_some_acids_multivariate[n,-c(1:3)])
    M <- cor_s_all[colnames (concentration_of_some_acids_multivariate[-c(1:3)]), colnames (concentration_of_some_acids_multivariate[-c(1:3)])]

    chi_concentration_of_some_acids_multivariate [n] <- t(Z) %*% solve(M) %*% Z

    }


p_value_concentration_of_some_acids_multivariate <- matrix(ncol=1, nrow=nrow(concentration_of_some_acids_multivariate))
p_value_concentration_of_some_acids_multivariate <- pchisq(chi_concentration_of_some_acids_multivariate, df=13, low=F)
rownames(p_value_concentration_of_some_acids_multivariate) <- rownames(locus_table)

save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_concentration_of_some_acids.Rdata", list=c("p_value_concentration_of_some_acids_multivariate"))

## 8th trait (pH)

ph_multivariate <- cbind(locus_table[,c("SNP", "CHR", "POS")], s_all_clear[,c("PH24LL", "PH24ST")])

chi_ph_multivariate <- matrix(ncol=1, nrow=nrow(ph_multivariate))
rownames(chi_ph_multivariate) <- rownames(locus_table)

for (n  in (1:nrow(ph_multivariate))) {

    Z <- as.numeric(ph_multivariate[n,-c(1:3)])
    M <- cor_s_all[colnames (ph_multivariate[-c(1:3)]), colnames (ph_multivariate[-c(1:3)])]

    chi_ph_multivariate [n] <- t(Z) %*% solve(M) %*% Z

    }


p_value_ph_multivariate <- matrix(ncol=1, nrow=nrow(ph_multivariate))
p_value_ph_multivariate <- pchisq(chi_ph_multivariate, df=13, low=F)
rownames(p_value_ph_multivariate) <- rownames(locus_table)

save(file="/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_ph.Rdata", list=c("p_value_ph_multivariate"))
