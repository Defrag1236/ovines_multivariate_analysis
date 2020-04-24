### make plot ###

# path to function
  source("/home/common/projects/ovine_selection/ovines_multivariate_analysis/zlobin_src/ovines_multivariate_analysis/05_make_regional_association_plots/05a_function_for_plot.R")

# load data 

new_snps <- read.table("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_snps_to_regional_association_plot.txt", head=F, stringsAsFactors=F)

for (n in (1:nrow(new_snps))) {
  # interesting (best) SNP
  snp       <- new_snps[n,1]
  # title of plot
  locusname <- paste("SNP", snp,"for multi_weight", sep=" ")
  # locus file (e.g. output of "building_locus_file.r")
  locus     <- read.table(paste("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot/locus_table_", snp, ".txt", sep=""), header=T, stringsAsFactors=F)
  # known genes path
	known_genes_path <- "na"	
  # path for plot
  pfad_ausgabe <- paste("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/plots/assocplot_LDregions_", snp,".pdf", sep="")
  # range of y-axis (P-values) 
  minimaler_pwert <- min(locus[,6])
  range_graph     <- round(-(log10(minimaler_pwert))+1)
  # pvalue of interesting (best) SNP
  pval_snp    <- locus[locus[,5]==snp,6]
  # Build dependend on used HapMap version for SNAP informations
  build_ncbi       <- "b36"
  # plattform used 
  plattform   <- "linux"

  # open plot
  pdf(pfad_ausgabe, width=6, height=6)
  # call function
  make.fancy.locus.plot(snp, locusname, locus, range_graph, pval_snp, build_ncbi, plattform, known_genes_path)
  # close plot
  dev.off()

}
