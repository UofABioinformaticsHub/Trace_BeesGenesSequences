library(dplyr)
library(Biostrings)

frag <- 

##Distance between 2 sequences
relativeDistance <- function(s1,s2){
  if (nchar(s1)!=nchar(s2)){
    pa <- pairwiseAlignment(s1,s2)
    s <- compareStrings(pa@pattern,pa@subject)
  } else {
    s <- compareStrings(s1,s2)
  }
  s <- strsplit(s,"")[[1]]
  sum(! (s  %in% c("A","T","G","C")))/length(s)
}

##Read the number of sequences in each sample
base <- paste0("AGRF_CAGRF15854_B6B68_MergedReads_xloose_readCounts/",frag)
f <- read.table(paste0(base,".txt"),stringsAsFactors = FALSE,header = FALSE) 
colnames(f) <- c("Count","Sequence")
nbReads <- read.table(paste0(base,".nb"))
nbReads$V2 <- gsub(paste0("(AGRF_CAGRF15854_B6B68_MergedReads_xloose_readCounts/|_B6B68_.+-.+-",frag,".txt)"),"",nbReads$V2)
colnames(nbReads) <- c("Count","Sample")
f[["Sample"]] <- rep(nbReads$Sample,nbReads$Count)
f <- as.data.frame(f)
f <- arrange(f,Sample,desc(Count))

##Cluster the sequences, select the major clusters for each sample and write them into .cl files  
d <- lapply(unique(f$Sample),function(sa){
    a <- filter(f,Sample==sa)
    a[["NbReads"]] <- sum(a$Count)
    
    i1 <- which(a$Count/sum(a$Count)>0.01)
    i2 <- which(a$Count>2)
    i <- intersect(i1,i2)
    if (length(i)==0) return(NULL) else ar <- a[i,]
    if (nrow(a)==1 && nrow(ar)==1){
      write.table(a,paste0("AGRF_CAGRF15854_B6B68_MergedReads_xloose_readClusters/",sa,"_",frag,".cl"),quote = FALSE,row.names = FALSE,col.names = TRUE)
      return(ar)
    }
    else {
      message(sa)
      check <- rep(0,nrow(a))
      for (i in 1:min(nrow(ar),nrow(a)-1)){
        if (check[i]==0){
          check[i] <- 1 
          for (j in (i+1):nrow(a)){
            dis <- relativeDistance(a$Sequence[i],a$Sequence[j])
            if (dis<=0.02){
              check[j] <- -1
              ar$Count[i] <- ar$Count[i]+a$Count[j]
            }
          }
        }
      }
      c <- ar[intersect(which(check==1),1:nrow(ar)),]
      c <- arrange(c,desc(Count))
      write.table(c,paste0("AGRF_CAGRF15854_B6B68_MergedReads_xloose_readClusters/",sa,"_",frag,".cl"),quote = FALSE,row.names = FALSE,col.names = TRUE)
      }
      return(c)
})
d <- do.call(rbind,d)

##Write the major clusters of each sample into a fasta file
dna <- DNAStringSet(d$Sequence)
names(dna) <- paste0(d$Sample,"_",d$Count)
writeXStringSet(dna,paste0("AGRF_CAGRF15854_B6B68_MergedReads_xloose_readClusters/",frag,".fasta"))

