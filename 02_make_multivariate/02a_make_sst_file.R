### make sst file ###

# install multiabel

pkgFile <- "/home/common/projects/Multivariate_analysis_IgG/soft/MultiABEL_1.1-6.tar.gz"

install.packages("svMisc") # dependencies for package

install.packages(pkgs=pkgFile, type="source", repos=NULL)

library(MultiABEL)

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("locus_table_from_Bolormaa.Rdata")

# make sst file

## make additional objectives

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/clear_trait_tables")

indep.snps <- locus_table[,1]
files.name <- c(1:61)


for (n in (6:61)) {

	x <- colnames(locus_table) [n]
	y <- paste (x, ".txt", sep="")
 
	files.name[n] <- y
	
	} 

files.name <- files.name[6:61]

## make sst

sst = load.summary(files=files.name,
	cor.pheno = NULL, indep.snps = as.character(indep.snps), est.var = TRUE,
	columnNames = c("MarkerName", "Allele1", "Allele0", "Freq1", "Effect", "StdErr", "N"), fixedN = NULL)

# save sst file 

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

save(sst, file="MV_sst_file.Rdata")