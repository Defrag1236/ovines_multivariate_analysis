### make clumping for the mv results steprens formula ###

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("locus_table_from_Bolormaa.Rdata")

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula")

load("multi_weight.Rdata")
load("multi_fat.Rdata")
load("multi_muscle.Rdata")
load("multi_concentration_of_some_substancies.Rdata")
load("multi_retail_colour.Rdata") 
load("multi_fresh_colour.Rdata")
load("multi_concentration_of_some_acids.Rdata")
load("multi_ph.Rdata")

library(data.table)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

shlop_all <- fread ("univariate_clumping.txt", head=T, stringsAsFactors=F, data.table=F)

# clumping function

function_for_shlop_28_12_2017=function(locus_table,p_value="P",pos="POS",snp="SNP",
                                       delta=2.5e5,chr="CHR",thr=5e-8,trait=NULL){
 #locus_table=bt
    locus_table[,p_value]=as.numeric(locus_table[,p_value])
    
    if (!is.null(trait)){
        traits="traits"
        locus_table=cbind(locus_table,traits=locus_table[,trait])
        locus_table[,traits]=as.character(locus_table[,traits])
    }
    
    out=locus_table[0,]
    
    locus_table=locus_table[locus_table[,p_value]<=thr,]
    
    i=1
    if (nrow(locus_table)>0){
        locus_table[,pos]=as.numeric(locus_table[,pos])
        locus_table[,p_value]=as.numeric(locus_table[,p_value])
        Zx <-locus_table
        Zx=Zx[order(Zx[,p_value]),]
        #n_traits=1
        #Zx=cbind(Zx,n_traits)
        i=1
        while (nrow(Zx)>0){       
            ind=which((abs(Zx[i,pos]-Zx[,pos])<=delta)&(Zx[i,chr]==Zx[,chr]))
                      
            if (!is.null(trait)){
                Zx[i,traits]=paste(unique(Zx[ind,trait]),collapse = ";")
            }
            
            out=rbind(out,Zx[i,])
            Zx=Zx[-ind,]            
        }
        rownames(out)=as.character(out[,snp])
    }
    
    if (!is.null(trait)){
        j=1
        out=cbind(out,Ntraits=1)
        out[,"Ntraits"]=as.numeric(out[,"Ntraits"])
        for (j in 1:nrow(out)){
            trs=unique(unlist(strsplit(out[j,traits],split = ";")))
            out[j,traits]=paste(trs,collapse = ";")
            out[j,"Ntraits"]=length(trs)
        }
    }
    
 return(out)
}

# do clumping for mv 

mv_stephens_list_p_value <- list(p_value_weight_multivariate, p_value_fat_multivariate, p_value_muscle_multivariate, p_value_concentration_of_some_substancies_multivariate, 
				p_value_retail_colour_multivariate, p_value_fresh_colour_multivariate, p_value_concentration_of_some_acids_multivariate, p_value_ph_multivariate)

mv_stephens_names_for_save <- c("weight_multivariate_stephens.txt", "fat_multivariate_stephens.txt", "muscle_multivariate_stephens.txt", "concentration_of_some_substancies_multivariate_stephens.txt", 
				"retail_colour_multivariate_stephens.txt", "fresh_colour_multivariate_stephens.txt", "concentration_of_some_acids_multivariate_stephens.txt", "p_value_ph_multivariate.txt")

names_mv <- c("multi_weight", "multi_fat", "multi_muscle",
    "multi_concentration_of_some_substancies", "multi_retail_colour", "multi_fresh_colour", "multi_concentration_of_some_acids", "multi_ph")

shlop_list_stephens_mv <- list()

m37 <- locus_table[,1:3]
colnames(m37) <- c ("SNP", "Chr", "Pos")

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_stephens_clumping")

for (n in (1:8)) {

	y <- mv_stephens_list_p_value[[n]]
	x <- names_mv[n]

    cycle_chr_pos <- cbind(m37, y, x)
    colnames(cycle_chr_pos) <- c("SNP", "CHR", "POS", "P", "trait")
    
    clumped <- function_for_shlop_28_12_2017(locus_table=cycle_chr_pos,p_value="P",pos="POS",snp="SNP",delta=5e5, chr="CHR",thr=0.05/nrow(locus_table), trait="trait")

  

    shlop_list_stephens_mv[[n]] <- clumped

    colnames(shlop_list_stephens_mv[[n]]) <- c("SNP", "CHR", "POS", "P", "trait", "traits", "Ntrait")

    fwrite(file=mv_stephens_names_for_save[n],x=clumped,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

	}

# do clumping for all mv stephens

shlop_all_stephens_mv <- do.call (rbind, shlop_list_stephens_mv)

shlop_result_all_stephens_mv <- function_for_shlop_28_12_2017(shlop_all_stephens_mv[,1:5], p_value="P", pos="POS",snp="SNP",
                                       delta=2.5e5,chr="CHR",thr=0.05/nrow(locus_table),trait="trait")

fwrite("clumping_for_all_mv_stephens.txt", x=shlop_result_all_stephens_mv,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

# do common clumping for Bolormaa and mv_stephens_Bolormaa

shlop_general_stephens <- rbind(shlop_all[,1:5], shlop_all_stephens_mv[,1:5])

shlop_result_general_stephens <- function_for_shlop_28_12_2017(shlop_general_stephens, p_value="P", pos="POS",snp="SNP",
                                       delta=2.5e5,chr="CHR",thr=0.05/nrow(locus_table),trait="trait")

fwrite("clumping_for_all_mv_stephens_and_uv.txt", x=shlop_result_general_stephens,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
