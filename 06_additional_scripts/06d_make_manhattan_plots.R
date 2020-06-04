### make manhattan plots ###

# load data

library(data.table)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

load("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_weight.Rdata")
load("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_fat.Rdata")
load("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula/multi_muscle.Rdata")

load("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata/locus_table_from_Bolormaa.Rdata")


# make tables for plot

chr_pos <- paste(snp_info$chromosome, snp_info$position, sep=":")

info_for_plot <- cbind(chr_pos, snp_info[,c(2,4:5)])
info_for_plot <- info_for_plot[!duplicated(info_for_plot$chr_pos),]

info_for_plot <- info_for_plot[which(info_for_plot[,1] %in% rownames(p_value_weight_multivariate)),]

rownames(info_for_plot) <- info_for_plot$chr_pos

p_value_weight_multivariate <- p_value_weight_multivariate[which(rownames(p_value_weight_multivariate) %in% rownames(info_for_plot)),]
p_value_fat_multivariate <- p_value_fat_multivariate[which(rownames(p_value_fat_multivariate) %in% rownames(info_for_plot)),]
p_value_muscle_multivariate <- p_value_muscle_multivariate[which(rownames(p_value_muscle_multivariate) %in% rownames(info_for_plot)),]

info_for_plot <- info_for_plot[order(rownames(info_for_plot), names(p_value_weight_multivariate)),]

for_plot_weight <- cbind(info_for_plot [,2:4], p_value_weight_multivariate)
colnames(for_plot_weight) <- c("SNP", "CHR", "BP", "P")

for_plot_fat <- cbind(info_for_plot [,2:4], p_value_fat_multivariate)
colnames(for_plot_fat) <- c("SNP", "CHR", "BP", "P")

for_plot_muscle <- cbind(info_for_plot [,2:4], p_value_muscle_multivariate)
colnames(for_plot_muscle) <- c("SNP", "CHR", "BP", "P")

# make chr column numeric

for_plot_weight$CHR <- as.numeric(for_plot_weight$CHR)
for_plot_fat$CHR <- as.numeric(for_plot_fat$CHR)
for_plot_muscle$CHR <- as.numeric(for_plot_muscle$CHR)

# make manhattan plots

library(qqman)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/plots/manhattan")

	## weight

	png("manhattan_weight_multivariate.png", width=1920, height=1080)

	manhattan(for_plot_weight, main="Multi weight manhattan plot", ylim=c(0,20), genomewideline=-log10(0.05/(nrow(locus_table)*(48+8))), suggestiveline=F, cex=0.8, 
    cex.axis=0.9, col=c("blue4", "orange3"))

	dev.off()

	## fat

	png("manhattan_fat_multivariate.png", width=1920, height=1080)

	manhattan(for_plot_fat, main="Multi fat manhattan plot", ylim=c(0,20), genomewideline=-log10(0.05/(nrow(locus_table)*(48+8))), suggestiveline=F, cex=0.8, 
    cex.axis=0.9, col=c("blue4", "orange3"))

	dev.off()
	
	## muscle

	png("manhattan_muscle_multivariate.png", width=1920, height=1080)

	manhattan(for_plot_muscle, main="Multi muscle manhattan plot", ylim=c(0,20), genomewideline=-log10(0.05/(nrow(locus_table)*(48+8))), suggestiveline=F, cex=0.8, 
    cex.axis=0.9, col=c("blue4", "orange3"))

	dev.off()


# make qq plots

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/plots/qqplots")

	## weight

	png("qq_multi_weigth.png", width=1920, height=1080)

	qq(for_plot_weight$P, main="Q-Q plot of multi weight p-values")

	dev.off()

	# fat

	png("qq_multi_fat.png", width=1920, height=1080)

	qq(for_plot_fat$P, main="Q-Q plot of multi fat p-values")

	dev.off()

	# muscle

	png("qq_multi_muscle.png", width=1920, height=1080)

	qq(for_plot_muscle$P, main="Q-Q plot of multi muscle p-values")

	dev.off()


