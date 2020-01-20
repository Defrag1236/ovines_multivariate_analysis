### make univariate clumping ###

# load data

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("locus_table_from_Bolormaa.Rdata")

# make univariate clumping 

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

shlop_list <- list()

for (n in(6:61)) {

	x <- colnames(locus_table[n])

	shlop_result <-  function_for_shlop_28_12_2017(locus_table[,c(1:3, n)], p_value=x, pos="POS",snp="SNP",
                                       delta=2.5e5,chr="CHR",thr=0.05/nrow(locus_table),trait=NULL)

	trait  <- rep (x, times=nrow(shlop_result))

	shlop_list[[n]] <- cbind(shlop_result, trait)

	colnames(shlop_list[[n]]) <- c("SNP", "CHR", "POS", "P", "trait")
	
	}


# do common clumping for univariate clumping

shlop_all <- do.call (rbind, shlop_list)

shlop_result_all <- function_for_shlop_28_12_2017(shlop_all, p_value="P", pos="POS",snp="SNP",
                                       delta=2.5e5,chr="CHR",thr=0.05/nrow(locus_table),trait="trait")


# save clumping result

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

write.table(shlop_result_all, "univariate_clumping.txt", row.names=F, col.names=T, quote=F)