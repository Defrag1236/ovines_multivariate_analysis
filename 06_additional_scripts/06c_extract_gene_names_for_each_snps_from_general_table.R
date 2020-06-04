### extract gene names for each snps from general table ###

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

snps <- read.table("general_table_for_paper.txt", head=T, stringsAsFactors=F)

gene_annot <- read.table("/home/common/projects/ovine_selection/Data/oar_3.1_gene_info/oar_3.1_gene_info.txt", head=T, stringsAsFactors=F)

library(diplyr)

# extract gene names

gene_names <- list()

for (n in (1:nrow(snps))) {

	x <- subset(gene_annot, gene_annot$CHR==snps$CHR[n], select=1:5)
	x <- subset(x, x$START>(snps$POS[n]-500000), select=1:5)
	x <- subset(x, x$START<(snps$POS[n]+500000), select=1:5)

	y <- subset(gene_annot, gene_annot$CHR==snps$CHR[n], select=1:5)
	y <- subset(y, y$STOP>(snps$POS[n]-500000), select=1:5)
	y <- subset(y, y$STOP<(snps$POS[n]+500000), select=1:5)

	z <- setdiff(y,x)

	x <- rbind(x,z)

	gene_names[[n]] <- as.vector(x$NAME)

}

gene_names <- gene_names[c(1:5,10)]

new_snps_with_gene_names <- snps[c(1:5,10), c(1:3)]
new_snps_with_gene_names$gene_names <- gene_names

new_snps_with_gene_names$gene_names <- vapply(new_snps_with_gene_names$gene_names, paste, collapse = ", ", character(1L))

# save result


write.table(new_snps_with_gene_names,"new_mv_snps_with_gene_name.txt", col.names=T, row.names=F, quote=F)
