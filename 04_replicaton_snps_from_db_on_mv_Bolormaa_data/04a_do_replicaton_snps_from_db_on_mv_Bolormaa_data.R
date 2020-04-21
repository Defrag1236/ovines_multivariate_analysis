### 04a do replicaton snps from db on mv Bolormaa data ###


# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("locus_table_from_Bolormaa.Rdata")

database_snps <- read.table("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/snps_from_database.txt", head=T, stringsAsFactors=F, sep="\t")
pc_snps <- read.table("/home/common/projects/ovine_selection/Data/PC_data_lead_SNP_from_Bolormaa.txt", head=F, stringsAsFactors=F, sep="\t")

# prepare data

locus_table <- locus_table[,c("PSWT", "HCWT", "LLWT", "TOP", "RND", "BONE", "LEGBONE", "LMY", "DRESSING", "IMF", "LLFAT", "CFATSCORE", "CCFAT", "HGRFAT", "CFAT5", "CEMW", "CEMD", "CEMA")]

pc_snps_clear <- paste(pc_snps[,3], pc_snps[,4], sep=":")

database_snps_clear <- paste(database_snps$CHR, database_snps$POS, sep=":")

database_snps_clear <- database_snps_clear[!(database_snps_clear %in% pc_snps_clear)]

# do replication 

snps_to_replication <- database_snps_clear[(database_snps_clear %in% rownames(locus_table))]

rep_snps <- list()

for (n in (1:length(snps_to_replication))) {

	x <- locus_table[which(snps_to_replication[n]==rownames(locus_table)),]

	for (i in (1:ncol(x))) {

		rep_snp <- matrix(ncol=3, nrow=1)

		colnames(rep_snp) <- c("SNP", "P", "trait")

		if (x[,i]<0.05/(18*length(database_snps_clear))) {

			
			
			rep_snp[1,1] <- rownames(x)
			rep_snp[1,2] <- min(x)
			rep_snp[1,3] <- colnames(x)[which.min(x)]
			
			rep_snps[[n]] <- rep_snp

			print("well done")
		}


	}

}

rep_snps <- do.call(rbind, rep_snps)


# save result

write.table(rep_snps, "/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/replicated_snps_from_db.txt", col.names=TRUE,row.names=FALSE,quote=FALSE)