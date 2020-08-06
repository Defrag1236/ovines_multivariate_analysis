### count explained dispression 

# load_data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("locus_table_from_Bolormaa.Rdata")

# prepare data

chr_pos <- strsplit(rownames(locus_table), ":")
chr_pos <- do.call(rbind, chr_pos)

chr_pos <- as.data.frame(chr_pos)
chr_pos[,1] <- as.numeric(as.character(chr_pos[,1]))
chr_pos[,2] <- as.numeric(as.character(chr_pos[,2]))

rownames(chr_pos) <- rownames(locus_table)

# count explained dispersion

	## 6 chr 31-42mb

	chr6 <- subset(chr_pos, chr_pos[,1]==6, select=1:2)
	chr6 <- subset(chr6, chr6[,2]>=31000000, select=1:2)
	chr6 <- subset(chr6, chr6[,2]<=42000000, select=1:2)

	chr6_p <- min(locus_table[which(rownames(locus_table) %in% rownames(chr6)),6:61])

	chr6_ed <- (qchisq(chr6_p, df=1, lower.tail=F)/10613)*100


	## 11 chr 24.5-28mb

	chr11 <- subset(chr_pos, chr_pos[,1]==11, select=1:2)
	chr11 <- subset(chr11, chr11[,2]>=24500000, select=1:2)
	chr11 <- subset(chr11, chr11[,2]<=28000000, select=1:2)

	chr11_p <- min(locus_table[which(rownames(locus_table) %in% rownames(chr11)),6:61])

	chr11_ed <- (qchisq(chr11_p, df=1, lower.tail=F)/10613)*100

	## 18 chr 62-67mb

	chr18 <- subset(chr_pos, chr_pos[,1]==18, select=1:2)
	chr18 <- subset(chr18, chr18[,2]>=62000000, select=1:2)
	chr18 <- subset(chr18, chr18[,2]<=67000000, select=1:2)

	chr18_p <- min(locus_table[which(rownames(locus_table) %in% rownames(chr18)),6:61])

	chr18_ed <- (qchisq(chr18_p, df=1, lower.tail=F)/10613)*100

# save result

all_ed <- data.frame("locus"=c("6_chr_31_42mb", "11chr_24.5_28mb", "18chr_62_67mb"), "ed_%"=rbind(chr6_ed, chr11_ed, chr18_ed)) 

write.table(all_ed, "/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/explainde_dispression_for_merged_loci.txt", col.names=T, row.names=F, quote=F)