### make list of snps for vep ###

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot")

snp1 <- read.table("locus_table_1:198273462.txt", head=T, stringsAsFactors=F)
snp2 <- read.table("locus_table_2:68158297.txt", head=T, stringsAsFactors=F)
snp3 <- read.table("locus_table_3:60513720.txt", head=T, stringsAsFactors=F)
snp4 <- read.table("locus_table_6:22846937.txt", head=T, stringsAsFactors=F)
snp5 <- read.table("locus_table_6:29898795.txt", head=T, stringsAsFactors=F)
snp6 <- read.table("locus_table_23:44492468.txt", head=T, stringsAsFactors=F)

library(data.table)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

# extract snps with r^2>=0.7

snp1 <- subset(snp1, snp1$RSquared>=0.7, select=5)
snp2 <- subset(snp2, snp2$RSquared>=0.7, select=5)
snp3 <- subset(snp3, snp3$RSquared>=0.7, select=5)
snp4 <- subset(snp4, snp4$RSquared>=0.7, select=5)
snp5 <- subset(snp5, snp5$RSquared>=0.7, select=5)
snp6 <- subset(snp6, snp6$RSquared>=0.7, select=5)

to_vep <- rbind(snp1, snp2, snp3, snp4, snp5, snp6) 

# extract rs_id

rs_id <- cbind(snp_info$rs, paste(snp_info$chromosome, snp_info$position, sep=":"))

to_vep_rs_id <- rs_id[which(rs_id[,2] %in% to_vep[,1]), 1]

# save results

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

write.table(to_vep, "snps_for_vep_r2_07.txt", col.names=F, row.names=F, quote=F)
write.table(to_vep_rs_id, "snps_for_vep_r2_07_rs_id.txt", col.names=F, row.names=F, quote=F)