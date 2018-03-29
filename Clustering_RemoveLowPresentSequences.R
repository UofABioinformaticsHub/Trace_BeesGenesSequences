source("utils.R")

#Clustering sequences for each sample in each gene fragment

frag <- 

base <- paste0(count_dir,frag)

#read the sequences file
f <- read.table(paste0(base,".txt"),stringsAsFactors = FALSE,header = FALSE)
colnames(f) <- c("Count","Sequence")

#read the counts file
nbReads <- read.table(paste0(base,".nb"))
nbReads$V2 <- gsub(paste0("(",count_dir,"|_",frag,".txt)"),"",nbReads$V2)
colnames(nbReads) <- c("Count","Species")

counts <- list()
for (s in nbReads$Species){
  nb <- read.table(paste0(count_dir,s,"_",frag,".nb"))
  names <- gsub(paste0("(",count_dir,"|_",frag,".txt)"),"",nb$V2)
  nb$V2 <- gsub("_.+$","",names)
  nb[["Sample"]] <- sapply(1:nrow(nb),function(i){gsub(paste0(nb$V2[i],"_"),"",names[i])})
  colnames(nb) <- c("Count","Species","Sample")
  counts[[s]] <- nb
}
counts <- do.call(rbind,counts)

f[["Species"]] <- rep(counts$Species,counts$Count)
f[["Sample"]] <- rep(counts$Sample,counts$Count)
f <- as.data.frame(f)
f <- arrange(f,Species,Sample,desc(Count))

allClusters <- lapply(unique(f$Sample),function(sa){
  sampleCluster <- filter(f,Sample==sa)
  sampleCluster[["NbReads"]] <- sum(sampleCluster$Count)
  
  #remove every sequence appears least than 3 times, or appears least than 1% of the sample
  i1 <- which(sampleCluster$Count/sum(sampleCluster$Count)>0.01)
  i2 <- which(sampleCluster$Count>2)
  i <- intersect(i1,i2)
  if (length(i)==0) return(NULL) else sampleKeptCluster <- sampleCluster[i,]
  
  #if there's only one sequence remaining, then write it into the clustered file, don't need to process the clustering
  if (nrow(sampleCluster)==1 && nrow(sampleKeptCluster)==1){
    write.table(sampleCluster,paste0(cluster_dir,sa,"_",frag,".cl"),quote = FALSE,row.names = FALSE,col.names = TRUE)
    return(sampleKeptCluster)
  }
  else {#process the clustering when we have more than 1 sequences remaining
    message(sa)
    check <- rep(0,nrow(sampleCluster)) #a vector to mark if a sequence is grouped to another sequence (check=-1), or kept as the main sequence of the cluster (check=1)
    for (i in 1:min(nrow(sampleKeptCluster),nrow(sampleCluster)-1)){
      if (check[i]==0){
        check[i] <- 1 
        for (j in (i+1):nrow(sampleCluster)){
          if (check[j]==0){
            dis <- relativeDistance(sampleCluster$Sequence[i],sampleCluster$Sequence[j])
            if (dis<=0.02){
              check[j] <- -1
              sampleKeptCluster$Count[i] <- sampleKeptCluster$Count[i]+sampleCluster$Count[j]
            }
          }
        }
      }
    }
    if (check[nrow(sampleCluster)]==0) check[nrow(sampleCluster)] <- 1
    sampleKeptCluster <- sampleKeptCluster[intersect(which(check==1),1:nrow(sampleKeptCluster)),] #keep only the main sequence of each cluster (with check=1)
    sampleKeptCluster <- arrange(sampleKeptCluster,desc(Count)) #arrange the sequences of the number of reads
    write.table(sampleKeptCluster,paste0(cluster_dir,sa,"_",frag,".cl"),quote = FALSE,row.names = FALSE,col.names = TRUE) #write the squences into cluster folder
  }
  return(sampleKeptCluster)
})
allClusters <- do.call(rbind,allClusters)
dna <- DNAStringSet(allClusters$Sequence)
names(dna) <- paste0(allClusters$Species,"_",allClusters$Sample,"_",allClusters$Count)
writeXStringSet(dna,paste0(cluster_dir,frag,".fasta"))

