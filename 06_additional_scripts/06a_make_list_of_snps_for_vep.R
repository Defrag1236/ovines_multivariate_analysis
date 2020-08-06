### make list of snps for vep ###

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot")

snp1 <- read.table("locus_table_1:27641947.txt", head=T, stringsAsFactors=F)
snp2 <- read.table("locus_table_1:40823201.txt", head=T, stringsAsFactors=F)
snp3 <- read.table("locus_table_1:51978516.txt", head=T, stringsAsFactors=F)
snp4 <- read.table("locus_table_1:59961282.txt", head=T, stringsAsFactors=F)
snp5 <- read.table("locus_table_1:198273462.txt", head=T, stringsAsFactors=F)
snp6 <- read.table("locus_table_2:68158297.txt", head=T, stringsAsFactors=F)
snp7 <- read.table("locus_table_2:69900354.txt", head=T, stringsAsFactors=F)
snp8 <- read.table("locus_table_3:60513720.txt", head=T, stringsAsFactors=F)
snp9 <- read.table("locus_table_3:153924034.txt", head=T, stringsAsFactors=F)
snp10 <- read.table("locus_table_3:219082890.txt", head=T, stringsAsFactors=F)
snp11 <- read.table("locus_table_4:85985834.txt", head=T, stringsAsFactors=F)
snp12 <- read.table("locus_table_5:59475661.txt", head=T, stringsAsFactors=F)
snp13 <- read.table("locus_table_6:68973672.txt", head=T, stringsAsFactors=F)
snp14 <- read.table("locus_table_9:79387844.txt", head=T, stringsAsFactors=F)
snp15 <- read.table("locus_table_21:7437422.txt", head=T, stringsAsFactors=F)
snp16 <- read.table("locus_table_23:44492468.txt", head=T, stringsAsFactors=F)

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
snp7 <- subset(snp7, snp7$RSquared>=0.7, select=5)
snp8 <- subset(snp8, snp8$RSquared>=0.7, select=5)
snp9 <- subset(snp9, snp9$RSquared>=0.7, select=5)
snp10 <- subset(snp10, snp10$RSquared>=0.7, select=5)
snp11 <- subset(snp11, snp11$RSquared>=0.7, select=5)
snp12 <- subset(snp12, snp12$RSquared>=0.7, select=5)
snp13 <- subset(snp13, snp13$RSquared>=0.7, select=5)
snp14 <- subset(snp14, snp14$RSquared>=0.7, select=5)
snp15 <- subset(snp15, snp15$RSquared>=0.7, select=5)
snp16 <- subset(snp16, snp16$RSquared>=0.7, select=5)


to_vep <- rbind(snp1, snp2, snp3, snp4, snp5, snp6, snp7, snp8, snp9, snp10, snp11, snp12, snp13, snp14, snp15, snp16) 

# extract rs_id

rs_id <- cbind(snp_info$rs, paste(snp_info$chromosome, snp_info$position, sep=":"))

to_vep_rs_id <- rs_id[which(rs_id[,2] %in% to_vep[,1]), 1]

# save results

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

write.table(to_vep, "snps_for_vep_r2_07.txt", col.names=F, row.names=F, quote=F)
write.table(to_vep_rs_id, "snps_for_vep_r2_07_rs_id.txt", col.names=F, row.names=F, quote=F)