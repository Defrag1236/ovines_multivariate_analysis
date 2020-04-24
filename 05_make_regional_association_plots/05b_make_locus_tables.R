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

#ld_snp <- ld_clear[which(ld_clear[,1] %in% snp[,5]),]

#r2_1 <- ld_snp[grepl(pattern="1:198273462", x=ld_snp[,1]),]
#r2_2 <- ld_snp[grepl(pattern="1:198273462", x=ld_snp[,2]),]

#snp[which(snp[,5] %in% r2_1[,2]),3] <- r2_1[,3]
#snp[which(snp[,5] %in% r2_2[,1]),3] <- r2_1[,3]



# save results

write.table(snp, "/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot/locus_table_1_198273462.txt", row.names=F, quote=F)


write.table(for_plot, "/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot/locus_table.txt", row.names=F, quote=F)


