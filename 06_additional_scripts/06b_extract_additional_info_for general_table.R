### extract additional info for general table ###

# load data

library(data.table)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)


library("readxl")

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/clear_excel_tables")

general_table <- read_excel("mv_clumping_results.xlsx")

general_table <- as.data.frame(general_table)

rownames(general_table) <- paste(general_table$CHR, general_table$POS, sep=":")

load("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_weight.Rdata")
load("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_fat.Rdata")
load("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_muscle.Rdata")

# extract info

## extract p-value for each mv trait 

p_value_fat <- p_value_fat_multivariate[which(rownames(p_value_fat_multivariate) %in% rownames(general_table)),]
p_value_fat <- p_value_fat[match(rownames(general_table), names(p_value_fat))]

p_value_muscle <- p_value_muscle_multivariate[which(rownames(p_value_muscle_multivariate) %in% rownames(general_table)),]
p_value_muscle <- p_value_muscle[match(rownames(general_table), names(p_value_muscle))]



## extract alleles

snp_info <- cbind(snp_info$rs, snp_info$Alleles_A_B_FORWARD, paste(snp_info$chromosome, snp_info$position, sep=":"))

snp_info <- snp_info[which(snp_info[,3] %in% rownames(general_table) ),]

alleles <- do.call(rbind, (strsplit(snp_info[,2], "/")))
colnames(alleles) <- c("ra", "ea")


# do general table

general_table <- cbind(general_table, snp_info[,1], p_value_fat, p_value_muscle, alleles)

colnames(general_table)[8] <- c("rs")

# save result

write.table(general_table, "/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/general_table_for_paper.txt", col.names=T, row.names=F, quote=F)