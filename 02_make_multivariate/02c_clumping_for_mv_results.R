### clumping for mv results ###

# load data

library(data.table)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("locus_table_from_Bolormaa.Rdata")

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

# cbind mv_results with chr+pos and do clumping

idnames_mv <- c("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_weight.Rdata",
    "/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_fat.Rdata",
    "/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_muscle.Rdata",
    "/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_concentration_of_some_substancies.Rdata",
    "/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_retail_colour.Rdata",
    "/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_fresh_colour.Rdata",
    "/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_concentration_of_some_acids.Rdata",
    "/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/multi_ph.Rdata"
     )


names_mv <- c("multi_weight", "multi_fat", "multi_muscle",
    "multi_concentration_of_some_substancies", "multi_retail_colour", "multi_fresh_colour", "multi_concentration_of_some_acids", "multi_ph")

names_for_save_mv <- c("multi_weight.txt", "multi_fat.txt", "multi_muscle.txt",
    "multi_concentration_of_some_substancies.txt", "multi_retail_colour.txt", "multi_fresh_colour.txt", "multi_concentration_of_some_acids.txt", "multi_ph.txt")


m37 <- locus_table[,1:3]
colnames(m37) <- c ("SNP", "Chr", "Pos")

shlop_list_mv <- list()

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_multiabel_clumping")

for (n in 1:8) {

    load(idnames_mv[n])
    eval(expr=parse(text=paste("MV_P=",names_mv[n])))
    x <- names_mv[n]
    eval(expr=parse(text=paste("rm(",names_mv[n],")")))

    y=MV_P$scan
   
    cycle_chr_pos <- y[(y[,1] %in% m37[,1]),]
    ind <- match(rownames(cycle_chr_pos),m37[,1])
    Chr <- m37[ind,2]
    Pos <- m37[ind,3]
    cycle_chr_pos <- cbind(cycle_chr_pos[,1:4], Chr, Pos, x)

    
    clumped <- function_for_shlop_28_12_2017(locus_table=cycle_chr_pos,p_value="p",pos="Pos",snp="marker",delta=5e5, chr="Chr",thr=0.05/nrow(locus_table), trait="x")

    clumped <- clumped[,-c(2:3,7,9)]

    clumped <- clumped[,c(1,3,4,2,5)]

    shlop_list_mv[[n]] <- clumped

    colnames(shlop_list_mv[[n]]) <- c("SNP", "CHR", "POS", "P", "trait")

    fwrite(file=names_for_save_mv[n],x=clumped,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
    }

# do clumping for all mv

shlop_all_mv <- do.call (rbind, shlop_list_mv)

shlop_result_all_mv <- function_for_shlop_28_12_2017(shlop_all_mv, p_value="P", pos="POS",snp="SNP",
                                       delta=2.5e5,chr="CHR",thr=0.05/nrow(locus_table),trait="trait")

fwrite("clumping_for_all_mv_multiabel.txt", x=shlop_result_all_mv,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

# do common clumping for Bolormaa and mv_Bolormaa

shlop_general <- rbind(shlop_all[,1:5], shlop_all_mv)


shlop_result_general <- function_for_shlop_28_12_2017(shlop_general, p_value="P", pos="POS",snp="SNP",
                                       delta=2.5e5,chr="CHR",0.05/nrow(locus_table),trait="trait")

fwrite("clumping_for_mv_multiabel_and_uv_Bolormaa_data.txt", x=shlop_result_general,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
