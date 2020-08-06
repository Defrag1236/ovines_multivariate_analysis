###count_r2_between_new_mv_snps and do graphs and Hierarchical clustering plot

#load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula")

load("multi_weight.Rdata")

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("supplement_from_Bolormaa.Rdata")

new_snps <- read.table("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/mv_snps_new_threshold.txt", head=F, stringsAsFactors=F)
gene_names <- read.table("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/mv_snps_gene_names.txt", head=F, stringsAsFactors=F)

#count r2

z_snp <- s_all[which(rownames(s_all) %in% new_snps[,1]),]

r_snp <- cor(t(z_snp), t(z_snp))

r2_snp <- r_snp^2 

# count r2_spearman and p-value

r_snp_s <- cor(t(z_snp), t(z_snp), method="spearman")

new_snps_z <- t(z_snp)

r2_snp_s <- r_snp_s^2 

cor_r2_p <- function(new_snps_z, ...) {
    mat <- as.matrix(new_snps_z)
    n <- ncol(new_snps_z)
    p_r2_s<- matrix(NA, n, n)
    diag(p_r2_s) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(new_snps_z[, i], new_snps_z[, j], method="spearman", ...)
            p_r2_s[i, j] <- p_r2_s[j, i] <- tmp$p.value
        }
    }
  colnames(p_r2_s) <- rownames(p_r2_s) <- colnames(new_snps_z)
  p_r2_s
}


r2_p_s <- cor_r2_p(t(z_snp))


# make plot

library(corrplot)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/plots")

png("r2_corrplot_mv_snps.png")

corrplot(r2_snp, p.mat=round(r2_snp,1), insig ="p-value", sig.level=-1, method="square")

dev.off()

# make plot with significance for spearman correlation

png("r2_corrplot_mv_snps_spearman_with_significance.png", width=1920, height=1080)

corrplot(r2_snp_s, insig ="pch", sig.level=0.05/((20*19)/2), method="square", p.mat=r2_p_s)

dev.off()

# make  Hierarchical clustering plot


r2_snp_s_with_gene_names <- r2_snp_s
colnames(r2_snp_s_with_gene_names) <- gene_names[,2]
rownames(r2_snp_s_with_gene_names) <- gene_names[,2]

png("Heatmap_mv_snps.png", width=1920, height=1080)

heatmap(r2_snp_s_with_gene_names, scale = "none")

dev.off()
