### make locus tables for regional association plot ### 


# load 

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula")

load("multi_weight.Rdata")

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("supplement_from_Bolormaa.Rdata")

new_snps <- read.table("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_snps_to_regional_association_plot.txt", head=F, stringsAsFactors=F)

#library(data.table)

#ld <- fread("/home/common/projects/ovine_selection/ovines_gwas_map/sheep_reference/sheep_reference.ld", head=T, stringsAsFactors=F, data.table=F)

# make chr_pos for locus table 

chr_pos <- strsplit(rownames(p_value_weight_multivariate), ":")
chr_pos <- do.call(rbind, chr_pos)
chr_pos<- apply(chr_pos, 2, as.numeric)

# update ld table for needed format

#ld_a <- paste(ld$CHR_A, ld$BP_A, sep=":")
#ld_b <- paste(ld$CHR_B, ld$BP_B, sep=":")
#ld_clear <- cbind.data.frame(ld_a, ld_b, ld$R2, stringsAsFactors=F)

# make common locus table

for_plot <- data.frame("Coordinate_HG18"=chr_pos[,2],"Chromosome"= chr_pos[,1],"RSquared"=0, "RecombinationRate"=0, "Proxy"=rownames(p_value_weight_multivariate),
						 "PVAL"=p_value_weight_multivariate, "TYPE"=rep("imputed",length(p_value_weight_multivariate)), stringsAsFactors=F)

# make locus table for each new snp


for (n in (1:nrow(new_snps))) {

	top <- new_snps[n,1]

	top_chr_pos <- strsplit(top, ":")
	top_chr_pos <- do.call(cbind, top_chr_pos)
	top_chr_pos <- apply(top_chr_pos, 2, as.numeric)

	snp <- subset(for_plot,  for_plot[,2]==top_chr_pos[1,1], select=1:7)
	snp <- subset(snp, snp[,1]<(top_chr_pos[2,1]+500000), select=1:7)
	snp <- subset(snp, snp[,1]>(top_chr_pos[2,1]-500000), select=1:7)

	z_snp <- s_all[which(rownames(s_all) %in% snp[,5]),]

	z_top <- z_snp[which(rownames(z_snp) %in% top),]

	r_snp <- cor(t(z_top), t(z_snp))

	r2_snp <- r_snp^2 

	snp[,3] <- t(r2_snp) 

	colnames(snp[,3]) <- c("RSquared")

	write.table(snp, paste("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot/locus_table_", top, ".txt", sep=""), row.names=F, quote=F)

}

# make locus table for 2 interesting regions

## 6 chr 31-42 mb
	
snp_6chr <- subset(for_plot,  for_plot[,2]==6, select=1:7)
snp_6chr <- subset(snp_6chr, snp_6chr[,1]<42000000, select=1:7)
snp_6chr <- subset(snp_6chr, snp_6chr[,1]>31000000, select=1:7)

top_6chr <- snp_6chr[which.min(snp_6chr$PVAL),5]

top_chr_pos <- strsplit(top, ":")
top_chr_pos <- do.call(cbind, top_chr_pos)
top_chr_pos <- apply(top_chr_pos, 2, as.numeric)


z_snp_6chr <- s_all[which(rownames(s_all) %in% snp_6chr[,5]),]

z_top_6chr <- z_snp_6chr[which(rownames(z_snp_6chr) %in% top_6chr),]

r_snp_6chr <- cor(t(z_top_6chr), t(z_snp_6chr))

r2_snp_6chr <- r_snp_6chr^2 

snp_6chr[,3] <- t(r2_snp_6chr) 

colnames(snp_6chr[,3]) <- c("RSquared")

write.table(snp, "/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot/locus_table_6chr_31_42mb.txt", row.names=F, quote=F)

## 6 chr 31-42 mb
	
snp_6chr <- subset(for_plot,  for_plot[,2]==6, select=1:7)
snp_6chr <- subset(snp_6chr, snp_6chr[,1]<42000000, select=1:7)
snp_6chr <- subset(snp_6chr, snp_6chr[,1]>31000000, select=1:7)

top_6chr <- snp_6chr[which.min(snp_6chr$PVAL),5]

top_chr_pos <- strsplit(top, ":")
top_chr_pos <- do.call(cbind, top_chr_pos)
top_chr_pos <- apply(top_chr_pos, 2, as.numeric)


z_snp_6chr <- s_all[which(rownames(s_all) %in% snp_6chr[,5]),]

z_top_6chr <- z_snp_6chr[which(rownames(z_snp_6chr) %in% top_6chr),]

r_snp_6chr <- cor(t(z_top_6chr), t(z_snp_6chr))

r2_snp_6chr <- r_snp_6chr^2 

snp_6chr[,3] <- t(r2_snp_6chr) 

colnames(snp_6chr[,3]) <- c("RSquared")

write.table(snp_6chr, "/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot/locus_table_6chr_31_42mb.txt", row.names=F, quote=F)

## 11 chr 24-28.5 mb
	
snp_11chr <- subset(for_plot,  for_plot[,2]==11, select=1:7)
snp_11chr <- subset(snp_11chr, snp_11chr[,1]<28500000, select=1:7)
snp_11chr <- subset(snp_11chr, snp_11chr[,1]>24000000, select=1:7)

top_11chr <- snp_11chr[which.min(snp_11chr$PVAL),5]

top_chr_pos <- strsplit(top, ":")
top_chr_pos <- do.call(cbind, top_chr_pos)
top_chr_pos <- apply(top_chr_pos, 2, as.numeric)


z_snp_11chr <- s_all[which(rownames(s_all) %in% snp_11chr[,5]),]

z_top_11chr <- z_snp_11chr[which(rownames(z_snp_11chr) %in% top_11chr),]

r_snp_11chr <- cor(t(z_top_11chr), t(z_snp_11chr))

r2_snp_11chr <- r_snp_11chr^2 

snp_11chr[,3] <- t(r2_snp_11chr) 

colnames(snp_11chr[,3]) <- c("RSquared")

write.table(snp_11chr, "/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot/locus_table_11chr_24_28,5mb.txt", row.names=F, quote=F)


#ld_snp <- ld_clear[which(ld_clear[,1] %in% snp[,5]),]

#r2_1 <- ld_snp[grepl(pattern="1:198273462", x=ld_snp[,1]),]
#r2_2 <- ld_snp[grepl(pattern="1:198273462", x=ld_snp[,2]),]

#snp[which(snp[,5] %in% r2_1[,2]),3] <- r2_1[,3]
#snp[which(snp[,5] %in% r2_2[,1]),3] <- r2_1[,3]



# save results

write.table(for_plot, "/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot/locus_table.txt", row.names=F, quote=F)


