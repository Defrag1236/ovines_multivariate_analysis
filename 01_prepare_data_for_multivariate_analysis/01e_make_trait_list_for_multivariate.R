### make trait list for multivariate ###

# load data 

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("locus_table_from_Bolormaa.Rdata")

load("b_se_locus_tables.Rdata")

# make b and se tables in order

b_locus_table <- b_locus_table[order(rownames(b_locus_table)),]
se_locus_table <- se_locus_table[order(rownames(se_locus_table)),]

# make clear trait tables

trait_list <- colnames(locus_table[6:61])

clear_tables_for_analysis <- list()


for (n in (1:56)) {

    x <- as.data.frame (matrix(ncol=7, nrow=nrow(locus_table)))
    rownames(x) <- rownames(locus_table)
    colnames(x) <- c("MarkerName", "Allele1", "Allele0", "Freq1", "Effect", "StdErr", "N")
    x [,1] <- locus_table$SNP
    x [,2] <- locus_table$A1   
    x [,3] <- locus_table$A0
    x [,4] <- 0.5
    x [,5] <- b_locus_table[rownames(b_locus_table) %in% rownames(x),n]
    x [,6] <- se_locus_table[rownames(se_locus_table) %in% rownames(x),n]
    x [,7] <- 10613

 	clear_tables_for_analysis[[n]] <- x    
    }


names(clear_tables_for_analysis) <- trait_list

# save clear tables

library(data.table)

for (n in (1:56)) {

	x <- clear_tables_for_analysis [[n]]
	y <- paste ("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/clear_trait_tables/", trait_list[n], ".txt", sep="")
 

	fwrite(x=x, file=y, row.names=F, quote=F,col.names=T)
}