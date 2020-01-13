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

# make correlation matrix and plot

cor_s_all <- cor(s_all)

library(corrplot)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/plots")

png("correlation_matrix_for_traits_from_Bolormaa.png")

corrplot(cor_s_all)

dev.off()

# save common data.frame as .Rdata

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

save(s_all, file="supplement_from_Bolormaa.Rdata")