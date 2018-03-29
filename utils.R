library(dplyr)
library(Biostrings)
library(magrittr)

#define some parameters
fragments <- c("M414F","M84F","M202F","M702F","M82F","G0605F","E393F","E577F","E783F","beewgFor","W158F")
genes <- c("CO15","CO13","EF","WNT")
fragNames <- c("M414F_M84F","M202F_M702F_M82F","G0605F_E393F_E577F_E783F","beewgFor_W158F")

#define some directories
count_dir <- "AGRF_CAGRF15854_B6B68_Xloose_readCounts/"
cluster_dir <- "AGRF_CAGRF15854_B6B68_Xloose_readClusters/"
merge_dir <- "AGRF_CAGRF15854_B6B68_Xloose_mergedReads/"
sampleSpecificRef_dir <- "AGRF_CAGRF15854_B6B68_Xloose_sampleSpecificRef/"
sampleSpecificRef_checkReplicates_dir <- "AGRF_CAGRF15854_B6B68_Xloose_sampleSpecificRef_checkReplicates/"

system(paste0("mkdir -p ",count_dir))
system(paste0("mkdir -p ",cluster_dir))
system(paste0("mkdir -p ",merge_dir))
system(paste0("mkdir -p ",sampleSpecificRef_dir))
system(paste0("mkdir -p ",sampleSpecificRef_checkReplicates_dir))

#Function calculates relative distance between 2 sequences s1,s2
relativeDistance <- function(s1,s2){
    pa <- pairwiseAlignment(s1,s2)
    s <- compareStrings(pa@pattern,pa@subject)
    s <- strsplit(s,"")[[1]]
    sum(! (s  %in% c("A","T","G","C")))/length(s)
}

#Function calculates distance between 2 equal length sequences s1,s2
distance <- function(s1,s2){
  s <- compareStrings(s1,s2)
  s <- strsplit(s,"")[[1]]
  sum(s=="?")
}

#function to check if s1 and s2 is overlapped or not
#minOver the minimum overlap region length
#maxOver the maximum overlap region length
#thershold the maximum number of allowed mismatches in the overlap region
isOverlap <- function(s1,s2,minOver,maxOver,threshold){
  overlap <- rep(FALSE,maxOver-minOver+1)
  rDistance <- rep(0,maxOver-minOver+1)
  for (i in minOver:maxOver){
    ends1 <- substr(s1,nchar(s1)-i+1,nchar(s1))
    heads2 <- substr(s2,1,i)
    rDistance[i-minOver+1] <- distance(ends1,heads2)
    if (rDistance[i-minOver+1]<=threshold) {
      overlap[i-minOver+1] <- TRUE
    }
  } 
  message(which.min(rDistance)+minOver-1," ",min(rDistance))
  n <- which(overlap==TRUE)
  if (length(n)>0){
    n <- n[length(n)]+minOver-1
    return(TRUE)
  }
  return(FALSE)
}

#function to combine s1, s2 if they overlap
#return empty string if s1, s2 are not overlap
merge <- function(s1,s2,minOver,maxOver,threshold){
  for (i in minOver:maxOver){
    ends1 <- substr(s1,nchar(s1)-i+1,nchar(s1))
    heads2 <- substr(s2,1,i)
    rDistance <- distance(ends1,heads2)
    if (rDistance<=threshold) {
      cs <- consensusString(DNAStringSet(c(ends1,heads2)),ambiguityMap=IUPAC_CODE_MAP)
      return(paste0(substr(s1,1,nchar(s1)-i),cs,substr(s2,i+1,nchar(s2))))
    }
  } 
  return("")
}

#function to cluster a set of sequences
#seqs: the set of sequences to cluster
#counts: the count of each sequence
#threshold: the relative similarity of 2 sequences to be clustered into one group
clustering <- function(seqs,counts,threshold){
  if (length(seqs)>1){
    check <- rep(0,length(seqs))
    for (i in 1:(length(seqs)-1)){
      if (check[i]==0){
        check[i] <- 1 
        for (j in (i+1):length(seqs)){
          dis <- relativeDistance(seqs[i],seqs[j])
          if (dis<=threshold){
            check[j] <- -1
            counts[[i]] <- counts[[i]]+counts[[j]]
          }
        }
      }
    }
    return(list("Sequence"=seqs[check==1],"Count"=counts[check==1]))
  }
  else return(list("Sequence"=seqs,"Count"=counts))
}

#function that merges the sequences coming from the same samples but 2 different neighbor gene fragment
#namein1 the base of the first gene fragment file name
#namein2 the base of the second gene fragment file name
#nameout the base name of the output file
#minOver the minimum overlap region length
#maxOver the maximum overlap region length
#threshold the allowed mismatches in the overlap region
#maxThreshold the maximum allowed mismatches in the overlap region
#Final a logical parameter indicating if we attain the the last gene fragment or not
merge2 <- function(namein1,namein2,nameout,minOver,maxOver,threshold,maxThreshold,Final){
  #read the first gene fragment file
  sequences1 <- readDNAStringSet(paste0(merge_dir,namein1,".fasta"))
  names1 <- names(sequences1)
  ss <- strsplit(names1,":")
  ss <- sapply(ss, function(s){s[length(s)]})
  ss <- strsplit(ss,"_")
  species1 <- sapply(ss,function(i){i[2]})
  samples1 <- sapply(ss,function(i){i[3]})
  nbReads1 <- as.integer(sapply(ss,function(i){i[4]}))
  
  #read the second gene fragment file
  sequences2 <- readDNAStringSet(paste0(merge_dir,namein2,".fasta"))
  names2 <- names(sequences2)
  ss <- strsplit(names2,":")
  ss <- sapply(ss, function(s){s[length(s)]})
  ss <- strsplit(ss,"_")
  species2 <- sapply(ss,function(i){i[2]})
  samples2 <- sapply(ss,function(i){i[3]})
  nbReads2 <- as.integer(sapply(ss,function(i){i[4]}))
  
  #refer all available samples
  allSamples <- unique(c(samples1,samples2))
  
  #merge sequences for each sample
  allMergedSeq <- lapply(allSamples,function(s){
    id1 <- which(samples1==s)
    id2 <- which(samples2==s)
    s1 <- sequences1[id1]
    s2 <- sequences2[id2]
    if (length(s1)>0 && length(s2)>0){#if sample s presents both 2 fragments
      isV <- matrix(rep(FALSE,length(s1)*length(s2)),nrow=length(s1)) # a logical matrix such that isV[i,j]=TRUE iff s1[i] and s2[j] overlaps
      th <- threshold # the initial number of mismatches allowed in the overlap region, can be increased upto maxThreshold if no overlap region is found
      while(sum(isV)==0 && th <= maxThreshold){
        for (i in 1:length(s1)){
          for (j in 1:length(s2)){
            if (nbReads1[id1[i]]==0 || nbReads2[id2[j]]==0) {
              isV[i,j] <- TRUE
            } else if (isOverlap(s1[i],s2[j],minOver,maxOver,th)){
              isV[i,j] <- TRUE
            }
          }
        }
        th <- th + 1
      }
      th <- th - 1
      
      if(sum(isV)==0){#if not every sequence of the first fragment that overlaps with any sequence of the second fragment, then take either one fragment and consider the other fragment is missing
        r1 <- DNAStringSet(s1)
        names(r1) <- paste0(names(s1),":",namein2,"_",species2[id2[1]],"_",samples2[id2[1]],"_0")
        r2 <- DNAStringSet(s2)
        names(r2) <- paste0(namein1,"_",species1[id1[1]],"_",samples1[id1[1]],"_0",":",names(s2))
        mergedSeq <- DNAStringSet(c(r1,r2))
      }else{#if there exists overlaped sequences, then merge the overlap sequences and removed the ones that do not overlap to any other
        k1 <- rowSums(isV)
        k2 <- colSums(isV)
        x1 <- which(k1>0)
        mergedSeq <- lapply (x1,function(l){
          lapply(which(isV[l,]>0),function(v){
            if (nbReads1[id1[l]]==0 || nbReads2[id2[v]]==0){
              r1 <- xscat(s1[l],s2[v])
            } else{
              r1 <- DNAStringSet(merge(s1[l],s2[v],minOver,maxOver,th))  
            }
            names(r1) <- paste0(names(s1[l]),":",names(s2[v]))
            r1
          }) %>% DNAStringSetList() %>% unlist()
        }) %>% DNAStringSetList() %>% unlist()
      }
    } else if (length(s2)==0) {#missing sample s in second fragment
      mergedSeq <- DNAStringSet(s1)
      names(mergedSeq) <- paste0(names(s1),":",namein2,"_",species1[id1],"_",samples1[id1],"_0")
    } else if (length(s1)==0) {#missing sample s in first fragment
      mergedSeq <- DNAStringSet(s2)
      names(mergedSeq) <- paste0(namein1,"_",species2[id2],"_",samples2[id2],"_0",":",names(s2))
    }
    if (Final){#attain to the last fragment of the gene
      #rename and extract number of reads
      nameSeq <- names(mergedSeq)
      nb <- list()
      for (k in 1:length(mergedSeq)){
        f1 <- strsplit(nameSeq[k],":")[[1]]
        ss1 <- strsplit(f1,"_")
        n1 <- as.integer(sapply(ss1,function(s){s[length(s)]}))
        nb[[k]] <- n1
        sp <- sapply(ss1,function(s){s[2]})[1]        
        names(mergedSeq)[k] <- paste0(sp,"_",s,"_",paste0(n1,collapse = "_"))  
      }
      nbFrag <- length(nb[[1]])      
      
      #clustering again the merged sequences
      clusterMergedSeq <- clustering(mergedSeq,nb,0.02)
      mergedSeq <- clusterMergedSeq$Sequence
      nb <- clusterMergedSeq$Count
      
      #take the maximal merged sequence (in term of number of reads) 
      toRemove <- sapply(1:length(nb),function(i){
        any(sapply((1:length(nb))[-i],function(j){
          all(nb[[i]] <= nb[[j]]) || (nbFrag==3 && nb[[i]][3]==nb[[j]][3] && nb[[i]][2]==0)
        }))
      })
      
      #remove the low coverage sequences which does not overlap with any of its neighbor 
      toKeep <- rep(TRUE,length(mergedSeq))  
      for (k in 1:length(mergedSeq)){
        if (nb[[k]][1]<10 && nb[[k]][2]==0) nb[[k]][1] <- 0
        if (nb[[k]][length(nb[[k]])]<10 && nb[[k]][length(nb[[k]])-1]==0) nb[[k]][length(nb[[k]])] <- 0
        if (length(nb[[k]])>=3){
          for (j in 2:(length(nb[[k]])-1)){
            if (nb[[k]][j-1]==0 && nb[[k]][j]<10 && nb[[k]][j+1]==0) nb[[k]][j] <- 0
          }
        }
        if (all(n1==0)) toKeep[k] <- FALSE
      } 
      toKeep <- toKeep*(!toRemove)
      if (sum(toKeep)>0) mergedSeq <- mergedSeq[toKeep==1] else mergedSeq <- NULL
    }
    mergedSeq
  }) %>% DNAStringSetList() %>% unlist()
  writeXStringSet(allMergedSeq,paste0(merge_dir,nameout,".fasta"))
}


#function that merges the sequences coming from the same samples but 2 different neighbor gene fragment and check the consistency of replicates to infer appropriate reference
#toCheckSpecies the species to check
#namein1 the base of the first gene fragment file name
#namein2 the base of the second gene fragment file name
#nameout the base name of the output file
#minOver the minimum overlap region length
#maxOver the maximum overlap region length
#threshold the allowed mismatches in the overlap region
#maxThreshold the maximum allowed mismatches in the overlap region
#Final a logical parameter indicating if we attain the the last gene fragment or not
merge3 <- function(toCheckSpecies,namein1,namein2,nameout,minOver,maxOver,threshold,maxThreshold,Final){
  #load the first fragment sequences
  sequences1 <- readDNAStringSet(paste0(merge_dir,namein1,"_",toCheckSpecies,".fasta"))
  names1 <- names(sequences1)
  ss <- strsplit(names1,":")
  ss <- sapply(ss, function(s){s[length(s)]})
  ss <- strsplit(ss,"_")
  samples1 <- sapply(ss,function(i){i[3]})
  nbReads1 <- as.integer(sapply(ss,function(i){i[4]}))
  
  #load the second fragment sequences
  sequences2 <- readDNAStringSet(paste0(merge_dir,namein2,"_",toCheckSpecies,".fasta"))
  names2 <- names(sequences2)
  ss <- strsplit(names2,":")
  ss <- sapply(ss, function(s){s[length(s)]})
  ss <- strsplit(ss,"_")
  samples2 <- sapply(ss,function(i){i[3]})
  nbReads2 <- as.integer(sapply(ss,function(i){i[4]}))
  
  #infer all available samples
  allSamples <- unique(c(samples1,samples2))
  #the samples belong to the species to Check
  toCheckSamples <- unique(c(samples1,samples2))
  
  allMergedSeq <- lapply(toCheckSamples,function(s){
    id1 <- which(samples1==s)
    id2 <- which(samples2==s)
    s1 <- sequences1[id1]
    s2 <- sequences2[id2]
    if (length(s1)>0 && length(s2)>0){
      isV <- matrix(rep(FALSE,length(s1)*length(s2)),nrow=length(s1))
      th <- threshold
      while(sum(isV)==0 && th <= maxThreshold){
        for (i in 1:length(s1)){
          for (j in 1:length(s2)){
            if (nbReads1[id1[i]]==0 || nbReads2[id2[j]]==0) {
              isV[i,j] <- TRUE
            } else if (isOverlap(s1[i],s2[j],minOver,maxOver,th)){
              isV[i,j] <- TRUE
            }
          }
        }
        th <- th + 1
      }
      th <- th - 1
      k1 <- rowSums(isV)
      k2 <- colSums(isV)
      if(sum(isV)==0){
        r1 <- DNAStringSet(s1)
        names(r1) <- paste0(names(s1),":",namein2,"_",toCheckSpecies,"_",samples2[id2[1]],"_0")
        r2 <- DNAStringSet(s2)
        names(r2) <- paste0(namein1,"_",toCheckSpecies,"_",samples1[id1[1]],"_0",":",names(s2))
        mergedSeq <- DNAStringSet(c(r1,r2))
      }else{
        x1 <- which(k1>0)
        mergedSeq <- lapply (x1,function(l){
          lapply(which(isV[l,]>0),function(v){
            if (nbReads1[id1[l]]==0 || nbReads2[id2[v]]==0){
              r1 <- xscat(s1[l],s2[v])
            } else{
              r1 <- DNAStringSet(merge(s1[l],s2[v],minOver,maxOver,th))  
            }
            names(r1) <- paste0(names(s1[l]),":",names(s2[v]))
            r1
          }) %>% DNAStringSetList() %>% unlist()
        }) %>% DNAStringSetList() %>% unlist()
      }
    } else if (length(s2)==0) {
      mergedSeq <- DNAStringSet(s1)
      names(mergedSeq) <- paste0(names(s1),":",namein2,"_",toCheckSpecies,"_",samples1[id1],"_0")
    } else if (length(s1)==0) {
      mergedSeq <- DNAStringSet(s2)
      names(mergedSeq) <- paste0(namein1,"_",toCheckSpecies,"_",samples2[id2],"_0",":",names(s2))
    }
    if (Final){
      #rename the merged sequences
      nameSeq <- names(mergedSeq)
      idKeep <- rep(TRUE,length(mergedSeq))
      for (k in 1:length(mergedSeq)){
        f1 <- strsplit(nameSeq[k],":")[[1]]
        ss1 <- strsplit(f1,"_")
        n1 <- as.integer(sapply(ss1,function(s){s[length(s)]}))
        names(mergedSeq)[k] <- paste0(toCheckSpecies,"_",s,"_",paste0(n1,collapse = "_"))  
      }
    }
    mergedSeq
  }) %>% DNAStringSetList() %>% unlist()
  if (Final){
    ss <- strsplit(names(allMergedSeq),"_")
    samples <- sapply(ss,function(s){s[2]})
    nbR <- lapply(ss,function(s){as.integer(s[3:length(s)])}) 
    nbR <- do.call(rbind,nbR) %>% as.data.frame() 
    nbR <- mutate(nbR,TotalNbReads=rowSums(nbR))    
    
    #build a data frame contains information of each sequence to check 
    toCheck <- mutate(nbR,"Id"=1:nrow(nbR),"IdC"=1:nrow(nbR),"Sequence"=as.character(allMergedSeq),"Sample"=samples,"SampleInGroup"=samples,"Count"=TotalNbReads) %>% arrange(desc(Count))
    
    #cluster all merged sequences among all samples of the species
    check <- rep(0,nrow(toCheck))
    for (i in 1:(nrow(toCheck)-1)){
      if (check[i]==0){
        check[i] <- 1 
        for (j in (i+1):nrow(toCheck)){
          if (check[j]==0){
            dis <- relativeDistance(toCheck$Sequence[i],toCheck$Sequence[j])
            if (dis<=0.02){
              check[j] <- -1
              toCheck$Count[i] <- toCheck$Count[i]+toCheck$Count[j]
              toCheck$SampleInGroup[i] <- paste0(toCheck$SampleInGroup[i],",",toCheck$SampleInGroup[j])
              toCheck$IdC[j] <- toCheck$Id[i]
            }  
          }
        }
      }
    }
    if (check[nrow(toCheck)]==0) check[nrow(toCheck)] <- 1
    
    clusterOfSpecies <- toCheck[check==1,] 
    
    #calculate number of samples that support each sequence, and arrange by this number, followed by the number of Count
    clusterOfSpecies <- mutate(clusterOfSpecies,"NbSupportedSamples" = sapply(strsplit(SampleInGroup,","),function(s){length(unique(s))})) %>% 
                        arrange(desc(NbSupportedSamples),desc(Count)) %>% 
                        mutate("Order" = 1:nrow(clusterOfSpecies)) %>% 
                        select(IdC,Count,NbSupportedSamples,Order)
    
    toCheck <- full_join(toCheck,clusterOfSpecies,by="IdC") %>% select(Sample,IdC,Id,NbSupportedSamples,Order,TotalNbReads)
    
    #select the reference sequence for each sample
    selectRef <- sapply(toCheckSamples,function(s){
      ds <- filter(toCheck,Sample==s)
      mo <- min(ds$Order)
      nameSeq <- max(ds$TotalNbReads[which(ds$Order==mo)])
      ds$Id[ds$TotalNbReads==nameSeq][1]
    })
    writeXStringSet(allMergedSeq[selectRef],paste0(merge_dir,nameout,"_",toCheckSpecies,".fasta")) 
  } else{
    writeXStringSet(allMergedSeq,paste0(merge_dir,nameout,"_",toCheckSpecies,".fasta"))  
  }
}
