### function for plot

make.fancy.locus.plot <- function(snp, locusname, locus, range_graph, pval_snp, build_ncbi = "b36", plattform = "linux", known_genes_path="na") {

  # Postion entsprechend verwendetem build bestimmen
  if(build_ncbi=="b35")
    {
        POS <- locus$Coordinate_HG17
    }else{
     if(build_ncbi=="b36")
       {
        POS <- locus$Coordinate_HG18}else{stop("Invalid selection for NCBI Build!!! Please choose b35 or b36!!!")
       }
    }
  
  # Position an locus anhängen
  temp  <- data.frame(locus,"POS"=POS)

  # nach Position sortieren
  locus <- temp[order(temp$POS),]
  
  # untersuchen ob meherer Chomosomen enthalten sind, wenn ja abbrechen
  chr   <- unique(substring(locus$Chromosome,1))
  if(length(chr)!=1)stop(paste(length(chr),"different chromosomes!!!",sep=" "))
  
  # Zeile mit dem zentralen SNP
  hit   <- locus[locus$Proxy==snp,]
  

  # size of the region
  min.pos    <- min(locus$POS) - 10000
  max.pos    <- max(locus$POS) + 10000
  size.pos   <- max.pos - min.pos
  center.pos <- min.pos + ( size.pos / 2 )
  center.100kb.pos <- round(center.pos / 100000) * 100000
  offset.100kb.pos <- round((size.pos/3) / 100000) * 100000
  cat(min.pos,max.pos,size.pos,"\n",
  center.100kb.pos,"\t",offset.100kb.pos,"\n")
  # recombination rate
  # is included in SNAP output
  keep.recomb <- locus$RecombinationRate
  #
  # genes in the region
  #
  if (build_ncbi=="b35")
     {
      # alternatives Einlesen der Geninformation abhängig vom Betriebssystem
        if(known_genes_path=="na" & plattform == "windows")
          {
             genelist <- read.table(paste("O:\\Projects\\Software\\LD_GWA_plots\\0_MIT_BROAD_known_genes\\known_genes_build35_050307_chr", chr, ".txt", sep=""), header=T)
          }
        if(known_genes_path=="na" & plattform == "linux")
          {
             genelist <- read.table(paste("/mnt/nas/gwa/Projects/Software/LD_GWA_plots/0_MIT_BROAD_known_genes/known_genes_build35_050307_chr", chr, ".txt", sep=""), header=T)
          }
		genelist <- read.table(known_genes_path, header=T)  
        # Liste der Gene im Locusbereich erstellen
        genes.in.locus <- subset(genelist, ( genelist$START > min.pos & genelist$START < max.pos ) | ( genelist$STOP > min.pos & genelist$STOP < max.pos) )
     }else{
        if(build_ncbi=="b36")
          {
            # alternatives Einlesen der Geninformation abhängig vom Betriebssystem
            if(known_genes_path=="na" & plattform == "windows")
             {
               genes <- read.delim(paste("O:\\Projects\\Software\\LD_GWA_plots\\1_UCSC_genes_081120\\RefSeq_genes\\knownGenes_UCSC_mar_2008_chr", chr, ".txt", sep=""), header=T)
             }
            if(known_genes_path=="na" & plattform == "linux")
             {
               genes <- read.table(paste("/home/common/projects/ovine_selection/Data/oar_3.1_gene_info/oar_3.1_gene_info_chr_", chr, ".txt", sep=""), header=T)
             }
            # Liste der Gene im Locusbereich erstellen
            #genes <- read.table(known_genes_path, header=T) 
			genelist       <- data.frame("START"=genes$START,"STOP"=genes$STOP,"SIZE"=abs(genes$START-genes$STOP),"STRAND"=rep("+", length(genes[,1])),"GENE"=genes$NAME)
            genes.in.locus <- subset(genelist, ( genelist$START > min.pos & genelist$START < max.pos ) | ( genelist$STOP > min.pos & genelist$STOP < max.pos) )
          }
    }
  
  #print(genes.in.locus)
  n.genes     <- dim(genes.in.locus)[[1]]
  n.genes.pos <- dim(genes.in.locus[genes.in.locus$STRAND == "+",])[[1]]
  cat(n.genes,": ",n.genes.pos,"+ STRAND","\n")
  
  #
  # range_graph of y-axis
  #
  # this dedicates 33% of the yaxis to the genes, labels, recomb rate
  offset <- ( range_graph * 4 / 3 ) - range_graph
  big.range_graph <- range_graph + offset 
  
  ystart.gene <- - offset
  ystart.recomb <- - offset + (big.range_graph / 8)
  cat(big.range_graph,"\n",offset,"\n")
  
  #
  # genotyped markers
  #
  markers.in.strong.ld   <- subset(locus, (row.names(locus) != snp & locus$RSquared >= 0.8 & locus$TYPE == "typed"))
  markers.in.moderate.ld <- subset(locus, (row.names(locus) != snp & locus$RSquared >= 0.5 & locus$RSquared < 0.8 & locus$TYPE == "typed"))
  markers.in.weak.ld     <- subset(locus, (row.names(locus) != snp & locus$RSquared >= 0.2 & locus$RSquared < 0.5 & locus$TYPE == "typed"))
  markers.not.in.ld      <- subset(locus, (row.names(locus) != snp & locus$RSquared<0.2 & locus$TYPE == "typed"))
  #markers.typed <- subset(locus, (row.names(locus) != snp & locus$TYPE == "typed"))
  
  #
  # imputed SNPs
  #
  imputed.in.strong.ld   <- subset(locus, (row.names(locus) != snp & locus$RSquared >= 0.8 & locus$TYPE == "imputed"))
  imputed.in.moderate.ld <- subset(locus, (row.names(locus) != snp & locus$RSquared >= 0.5 & locus$RSquared < 0.8 & locus$TYPE == "imputed"))
  imputed.in.weak.ld     <- subset(locus, (row.names(locus) != snp & locus$RSquared >= 0.2 & locus$RSquared < 0.5 & locus$TYPE == "imputed"))
  imputed.not.in.ld      <- subset(locus, (row.names(locus) != snp & locus$RSquared<0.2 & locus$TYPE == "imputed"))
  #markers.imputed <- subset(locus, (row.names(locus) != snp & locus$TYPE == "imputed"))
  
  
  par(mar=c(4,4,3,4))
  
  #
  # start plot with recombination rate (in background)
  #
  plot(locus$POS, ystart.recomb + ( ( keep.recomb / 60 ) * ( 6 * big.range_graph / 8 )), type="l", col="lightblue", lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range_graph), xlab="", ylab="", main=locusname, axes=F)
  
  
  #
  # axes, titles and legends
  #
  mtext(paste("Chromosome", chr, "position (kb)", sep=" "), side=1, line=2.5)
  axis(1, at=c(center.100kb.pos - offset.100kb.pos, center.100kb.pos, center.100kb.pos + offset.100kb.pos), labels=c((center.100kb.pos - offset.100kb.pos) / 1000, center.100kb.pos / 1000, (center.100kb.pos + offset.100kb.pos) / 1000), las=1) 
  
  axis(2, at=seq(0,range_graph,2), labels=seq(0,range_graph,2), las=1)
  mtext("Observed (-logP)", side=2, at=(range_graph/2), line=2)
  
  axis(4, at=c( ystart.recomb, ystart.recomb + (big.range_graph / 4), ystart.recomb + ( 2 * big.range_graph / 4), ystart.recomb + ( 3 * big.range_graph / 4 ) ), labels=c("0","20","40","60"), las=1)
  mtext("Recombination rate (cM/Mb)", side=4, at=(-offset+big.range_graph/2), line=2)
  
  
  box()
  lines(c(min.pos, max.pos), c(0,0), lty="dotted", lwd=1, col="black")
  
  # zentralen SNP kennzeichnen
  if ( -(log10(pval_snp)) < range_graph )
    {
  	points(hit$POS, -(log10(pval_snp)), pch=23, cex=2.5, bg="blue")
  	text(hit$POS, -(log10(pval_snp)), labels=c(paste("P=",round(pval_snp, 15),sep="")), pos=4, offset=2)
  	text(hit$POS, -(log10(hit$PVAL)), labels=snp, pos=2, offset=1)
    }else{
  	points(hit$POS, range_graph, pch=23, cex=2.5, bg="blue")
  	text(hit$POS, range_graph, labels=c(paste("P=",round(pval_snp, 15),sep="")), pos=4, offset=1)
  	text(hit$POS, range_graph, labels=snp, pos=2, offset=1)
    }


  # besten SNP einzeichnen
  pos_best  <- locus$POS[locus$PVAL==min(locus$PVAL, na.rm=T)& is.na(locus$PVAL)==F]
  pval_best <- min(locus$PVAL, na.rm=T)
  snp_best  <- locus$Proxy[locus$PVAL==min(locus$PVAL, na.rm=T)& is.na(locus$PVAL)==F]
  if ( -(log10(pval_best)) < range_graph ) {
  	points(pos_best, -(log10(pval_best)), pch=23, cex=2.5, bg="green")
  	text(pos_best, -(log10(pval_best)), labels=c(paste("P=",round(pval_best, 15),sep="")), pos=4, offset=2)
  	text(pos_best, -(log10(pval_best)), labels=snp_best, pos=2, offset=1)
  } else {
  	points(pos_best, range_graph, pch=23, cex=2.5, bg="green")
  	text(pos_best, range_graph, labels=c(paste("P=",round(pval_best, 15),sep="")), pos=4, offset=1)
  	text(pos_best, -(log10(pval_best)), labels=snp_best, pos=2, offset=1)
  }
  
  #
  # plot the genotyped markers
  #
  points(markers.not.in.ld$POS, -(log10(markers.not.in.ld$PVAL)), pch=23, cex=1.0, bg="white")
  points(markers.in.weak.ld$POS, -(log10(markers.in.weak.ld$PVAL)), pch=23, cex=1.25, bg="yellow")
  points(markers.in.moderate.ld$POS, -(log10(markers.in.moderate.ld$PVAL)), pch=23, cex=1.25, bg="orange")
  points(markers.in.strong.ld$POS, -(log10(markers.in.strong.ld$PVAL)), pch=23, cex=1.25, bg="red")

  #
  # plot the imputed SNPs
  #
  points(imputed.not.in.ld$POS, -(log10(imputed.not.in.ld$PVAL)), pch=21, cex=1.0, bg="white")
  points(imputed.in.weak.ld$POS, -(log10(imputed.in.weak.ld$PVAL)), pch=21, cex=1.0, bg="grey")
  points(imputed.in.moderate.ld$POS, -(log10(imputed.in.moderate.ld$PVAL)), pch=21, cex=1.0, bg="orange")
  points(imputed.in.strong.ld$POS, -(log10(imputed.in.strong.ld$PVAL)), pch=21, cex=1.0, bg="red")

  #
  # plot the genes
  #---------------------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------------
  if (dim(genes.in.locus)[1]>0){
  # verändert am 18.02.09 von Janina Ried
  # Änderung: Gene anders gezeichnet, dazu
  #           - Bereiche identifizieren, in denen Gene überlappen
  #           - Zeichenbereich für die Gene wird in soviele Ebenen unterteilt wie es maximal überlappungen gibt
  #           - alle Gene werden abwechselnd auf diese Bereichen eingezeichnet
  
  # Gene sortieren erst nach der Stopposition, dann nach der Startposition
  genes.in.locus <- genes.in.locus[order(genes.in.locus$STOP),]
  genes.in.locus <- genes.in.locus[order(genes.in.locus$START),]
  # Doppelte Zeilen entfernen
  genes.in.locus <- genes.in.locus[duplicated(genes.in.locus)==F,]
  # unterschiedliche Splice-Varianten eines Gens werden zu einem Genbereich zusammengefasst
  # der Genbereich reicht vom Minimum der Startpositionen und das Maximum der Endpositionen 
  k <- 1
  klasse <- rep(-5, times=nrow(genes.in.locus))
  klasse[1] <- k
  if(nrow(genes.in.locus)>=2)
   {
      for(i in 2:nrow(genes.in.locus))
         {
           if(genes.in.locus$GENE[i]==genes.in.locus$GENE[i-1]){klasse[i] <- k}
           else{k <- k+1
                klasse[i] <-k}
         }
      genesections.in.locus <- matrix(-5, nrow=k,ncol=5)

      for(i in 1:k)
         {
           genesections.in.locus[i,1] <- min(genes.in.locus[klasse==i,1])
           genesections.in.locus[i,2] <- max(genes.in.locus[klasse==i,2])
           genesections.in.locus[i,3] <- as.integer(genesections.in.locus[i,2])-as.integer(genesections.in.locus[i,1])
           genesections.in.locus[i,4] <- levels(as.factor(as.character(genes.in.locus[klasse==i,4])))
           genesections.in.locus[i,5] <- levels(as.factor(as.character(genes.in.locus[klasse==i,5])))
         }
      colnames(genesections.in.locus) <- colnames(genes.in.locus)
      genes.in.locus <- as.data.frame(genesections.in.locus)
      genes.in.locus[,1]<- as.numeric(as.character(genes.in.locus[,1]))
      genes.in.locus[,2]<- as.numeric(as.character(genes.in.locus[,2]))
      genes.in.locus[,3]<- as.numeric(as.character(genes.in.locus[,3]))
   }
  print(genes.in.locus)
  
  # wenn mehr sehr kurze Gene aufeinander folgen, als es Stufen gibt kann es zu Problemen bei den Namen kommen, deshalb
  # zusätzliche Stufen einführen bis alle Gene einer Stufenfolge eine mindestlänge erfüllen
  # dazu bekommen Gene eine Mindestlänge zugeteilt
  # berechnet die Länge des Namens
  lab.width <- strwidth(genes.in.locus$GENE, units="user", cex=0.4)
  STOP_hyp <- genes.in.locus$START+2*lab.width
  #STOP_hyp <- genes.in.locus$START+20000
  STOP_max <- rep(1,length(genes.in.locus[,1]))
  for(i in 1:length(genes.in.locus[,1]))
     {
       STOP_max[i] <- max(genes.in.locus$STOP[i], STOP_hyp[i]) # mindestens die Länge des Namens
     }


  if(nrow(genes.in.locus)>=2)
    {
      diff <- rep(0,nrow(genes.in.locus))       # Differenz der Starposition und der Stopposition des vorherigen (positiv -> Überlappung)
      # Ueberlappung zwischen zwei Genen
      for(i in 2:nrow(genes.in.locus))
         {
           diff[i]<- STOP_max[i-1]-genes.in.locus$START[i]   # wenn positiv, dann ueberlappen sich die gene
         }
      k<-1
      # Zaehler der Ueberlappungen
      zaehler_ueberlap <- rep(1,nrow(genes.in.locus)) # zählt tatsächlich Überlappungen
      for (i in 2:nrow(genes.in.locus))
          {
            if(diff[i]>0)
              {
                k <- k+1                    # je mehr ueberlappungen desto mehr stufen
                zaehler_ueberlap[i] <- k
              }else{
                k=1
                zaehler_ueberlap[i] <- k
              }
          }

      # notwendige Stufen bestimmen
      # d.h. gibt es überlappung zwischen den Genen einer Stufen folge
      # nicht direkt aufeinanderfolgende Gene werden auf Überlappung untersucht
      # evtl. mögliche reduzierung der Stufenanzahl
      max_k       <- max(zaehler_ueberlap)        # höchste Stufe, die vorkommt
      if(max_k > 2)
        {
           haeuf_k     <- table(zaehler_ueberlap)
           haeuf_max_k <- haeuf_k[max_k]               # anzahl höchste Stufe
           stellen <- rep(0,haeuf_max_k)
           i=0
           for(j in 1:length(zaehler_ueberlap))
              {
                if( zaehler_ueberlap[j]==max_k)
                  {
                    i=i+1
                    stellen[i]=j
                  }
              }

          # schleife über alle Vorkommen der höchsten Stufe
          max_k_neu <- rep(max_k,haeuf_max_k)
          for(j in 1:haeuf_max_k)
          {
            stelle <- stellen[j]
            for(l in (max_k-1):2)
               {
                 if(STOP_max[stelle-l]-genes.in.locus$START[stelle] < 0){max_k_neu[j] <- max_k_neu[j] -1}
               }

          }

          # höchsten 8 Stufen sonst zu eng
          anzahl_stufen <- min(max(max_k_neu),8)
        }else{
          anzahl_stufen <- min(max(zaehler_ueberlap),8)
        }

      # Vektor erstellen in dem die Stufen eingetragen sind, auf denen die entsprechenden Gene eingezeichnet werden sollen
      zaehler_stufen <- seq(1,anzahl_stufen)                             # vorkommende Stufen
      zaehler        <- rep(zaehler_stufen,length.out=length(zaehler_ueberlap))  # abwechselnd diese Stufen eintragen

     }else{
     # falls nur eine Stufe vorkommt
      anzahl_stufen <- 1
      zaehler       <- 1
     }
  
  yteil.gene <- (big.range_graph/8)/anzahl_stufen # anteiliger Zeichenbereich für jede Stufe


  # Gene als Pfeile einzeichnen
  for( i in 1:nrow(genes.in.locus) )
    {
        for(ebene in 1:anzahl_stufen)
           {
                  if(zaehler[i]==ebene)
                    {
          	       if( genes.in.locus[i,]$STRAND == "+" )
                         {
                          arrows(max(genes.in.locus[i,]$START, min.pos), ystart.gene + (ebene-1)*yteil.gene, min(genes.in.locus[i,]$STOP, max.pos), ystart.gene + (ebene-1)*yteil.gene, length=0.05, lwd=2, code=2, lty="solid", col="darkgreen")
    	                 }else{
    		          arrows(max(genes.in.locus[i,]$START, min.pos), ystart.gene + (ebene-1)*yteil.gene, min(genes.in.locus[i,]$STOP, max.pos), ystart.gene + (ebene-1)*yteil.gene, length=0.05, lwd=2, code=1, lty="solid", col="darkgreen")
                         }

    	               if( ! is.na(genes.in.locus[i,]$GENE) )
                         {
                          #text(ifelse(genes.in.locus[i,]$START<min.pos,min.pos,min(max.pos,genes.in.locus[i,]$START + (genes.in.locus[i,]$SIZE / 2))), ystart.gene+ big.range_graph/80 + (ebene-1)*yteil.gene, labels=genes.in.locus[i,]$GENE, cex=0.5)
    	                  #text(ifelse(genes.in.locus[i,]$START<min.pos,min.pos,min(max.pos,genes.in.locus[i,]$START )), ystart.gene+ big.range_graph/80 + (ebene-1)*yteil.gene, labels=genes.in.locus[i,]$GENE, adj=c(0,0.5) ,cex=0.4)
                          text(max(genes.in.locus[i,]$START, min.pos)+(1/2)*(min(genes.in.locus[i,]$STOP, max.pos)-max(genes.in.locus[i,]$START, min.pos)), ystart.gene+ big.range_graph/80 + (ebene-1)*yteil.gene, labels=genes.in.locus[i,]$GENE, adj=c(0.5,0.5) ,cex=0.4)

                         }
                    }
           }
    }
  
  }
  #---------------------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------------  

}
