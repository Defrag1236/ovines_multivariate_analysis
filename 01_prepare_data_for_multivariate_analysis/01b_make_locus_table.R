### make locus table with p_value and correct SNP names ###

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("supplement_from_Bolormaa.Rdata")

setwd("/home/common/projects/ovine_selection/old/20181112")

load("20181113_snp_info_GG_num.RData")

# count p_value from z_score

s_all_p_value <- matrix(nrow=nrow(s_all), ncol=ncol(s_all))
rownames(s_all_p_value) <- rownames(s_all)
colnames(s_all_p_value) <- colnames(s_all)

for (n in (1:ncol(s_all))) {

	s_all_p_value[,n] <- pchisq (s_all[,n]^2, df=1, low=F)

	}

s_all_p_value <- s_all_p_value[order(rownames(s_all_p_value)),]

# make clear locus table

library(dplyr)
library(tidyr)

## make SNP column 

snp_name_from_chip <- cbind(snp_info[,1], paste(snp_info$Chr, snp_info$Position, sep=":"))
snp_name_from_chip <- snp_name_from_chip[snp_name_from_chip[,2] %in% rownames(s_all_p_value),]

snp_name_from_chip <- snp_name_from_chip[order(snp_name_from_chip[,2]),]

## separate and make chr and pos columns 

chr_pos <- as.data.frame(rownames (s_all_p_value))
chr_pos[,1] <- gsub ( ":", " ", x=chr_pos[,1])


chr_pos_clear <- chr_pos %>% separate (1, c("CHR", "POS"), sep="[\\s]+")
rownames(chr_pos_clear) <- rownames(s_all_p_value)

## compare and make the same dimensions for all columns

chr_pos_clear <- chr_pos_clear[rownames(chr_pos_clear) %in% snp_name_from_chip[,2],]
s_all_p_value <- s_all_p_value[rownames(s_all_p_value) %in% snp_name_from_chip[,2],]
snp_name_from_chip <- snp_name_from_chip[!duplicated(snp_name_from_chip[,2]),]


## rbind all in locus table

locus_table <- cbind(snp_name_from_chip[,1], chr_pos_clear, s_all_p_value)
colnames(locus_table) <- c("SNP", "CHR", "POS", colnames(s_all))


# save clear locus table as .Rdata

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

save(locus_table, file="locus_table_from_Bolormaa.Rdata")