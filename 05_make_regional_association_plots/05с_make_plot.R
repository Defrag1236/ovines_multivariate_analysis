### make plot ###

# path to function
source("/home/common/projects/ovine_selection/ovines_multivariate_analysis/zlobin_src/ovines_multivariate_analysis/05_make_regional_association_plots/05a_function_for_plot.R")

# load data 

new_snps <- read.table("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_snps_to_regional_association_plot.txt", head=F, stringsAsFactors=F)
chr6 <- read.table("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot/locus_table_6chr_31_42mb.txt", head=T, stringsAsFactors=F)
chr11 <-  read.table("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot/locus_table_11chr_24_28,5mb.txt", head=T, stringsAsFactors=F)
chr18 <-  read.table("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/for_regional_association_plot/locus_table_18chr_62_67mb.txt", head=T, stringsAsFactors=F)

# make plots for new snps associated with multi_weight

for (n in (c(1:8, 10, 13:20))) {
  # interesting (best) SNP
  snp       <- new_snps[n,1]
  # title of plot
  locusname <- paste("SNP", snp,"for MMass", sep=" ")
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


# make plots for new snps associated with multi_fat

for (n in (c(9, 11:12))) {
  # interesting (best) SNP
  snp       <- new_snps[n,1]
  # title of plot
  locusname <- paste("SNP", snp,"for MFat", sep=" ")
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

# make plots fo interesring regions

## 6chr 

  snp       <- chr6[which.min(chr6$PVAL),5]
  # title of plot
  locusname <- paste("6chr_31_42mb", sep="")
  # locus file (e.g. output of "building_locus_file.r")
  locus     <- chr6
  # known genes path
  known_genes_path <- "na"  
  # path for plot
  pfad_ausgabe <- paste("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/plots/assocplot_LDregions_6chr_31_42mb.pdf", sep="")
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

## 11chr 

  snp       <- chr11[which.min(chr11$PVAL),5]
  # title of plot
  locusname <- paste("11chr_24_28,5mb", sep="")
  # locus file (e.g. output of "building_locus_file.r")
  locus     <- chr11
  # known genes path
  known_genes_path <- "na"  
  # path for plot
  pfad_ausgabe <- paste("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/plots/assocplot_LDregions_11chr_24_28,5mb.pdf", sep="")
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

  ## 18 chr

    snp       <- chr18[which.min(chr18$PVAL),5]
  # title of plot
  locusname <- paste("18chr_62_67mb", sep="")
  # locus file (e.g. output of "building_locus_file.r")
  locus     <- chr18
  # known genes path
  known_genes_path <- "na"  
  # path for plot
  pfad_ausgabe <- paste("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/plots/assocplot_LDregions_18chr_62_67mb.pdf", sep="")
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
