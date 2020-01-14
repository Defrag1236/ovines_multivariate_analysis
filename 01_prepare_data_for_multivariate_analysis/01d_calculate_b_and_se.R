### calculate b and se ###

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("supplement_from_Bolormaa.Rdata")

load("locus_table_from_Bolormaa.Rdata")

# delete bad snps from z table

s_all <- s_all[rownames(s_all) %in% rownames(locus_table),]

# calculate se and b from z-score 

N <- 10613
p <- 0.5
q <- 0.5

z_locus_table <- s_all

se_locus_table <- apply (z_locus_table, 2, function(x) sqrt(1/(2*p*q*N*(1+x^2/N))))

b_locus_table <- matrix(ncol=ncol(z_locus_table), nrow=nrow(z_locus_table))

for (n in (1:nrow(z_locus_table))) {

	for (x in (1:ncol(z_locus_table))) {

		b_locus_table[n,x] <- se_locus_table[n,x]*z_locus_table[n,x] 

		}
	}

rownames(se_locus_table) <- rownames(s_all)
rownames(b_locus_table) <- rownames(s_all)


# save b and se tables

save (se_locus_table, b_locus_table, file="b_se_locus_tables.Rdata")