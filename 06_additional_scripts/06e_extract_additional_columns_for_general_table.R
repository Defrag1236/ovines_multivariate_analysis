### extract additional columns for general table ###

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("locus_table_from_Bolormaa.Rdata")

library("readxl")

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/clear_excel_tables")

general_table <- read_excel("mv_clumping_results.xlsx")

general_table <- as.data.frame(general_table)

rownames(general_table) <- paste(general_table$CHR, general_table$POS, sep=":")

# extract info

top_traits <- matrix(ncol=2, nrow=nrow(general_table))


for (n in 1:nrow(general_table)) {
	
	y <- rownames(general_table)[n]
	i <- locus_table[grepl(pattern=y, x=rownames(locus_table)),6:61]
	z <- which.min(i)
	top_traits[n,1] <- names(z)
	top_traits[n,2] <- i[,z]
}

# save result

write.table(top_traits, "/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/top_uv_traits_for_general_table.txt", col.names=F, row.names=F, quote=F) 
