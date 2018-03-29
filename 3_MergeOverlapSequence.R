source("utils.R")

#move the cluster files into the merged file and add the fragment name infront of each sequence name
for (frag in fragments){
  sequences <- readDNAStringSet(paste0(cluster_dir,frag,".fasta"))
  names(sequences) <- paste0(frag,"_",names(sequences))
  writeXStringSet(sequences,paste0(merge_dir,frag,".fasta"))
}

#merge neighbor gene fragments from the cluster files
merge2("M414F","M84F","M414F_M84F",7,20,2,3,TRUE)
merge2("M202F","M702F","M202F_M702F",60,70,2,10,FALSE)
merge2("M202F_M702F","M82F","M202F_M702F_M82F",20,30,2,5,TRUE)
merge2("G0605F","E393F","G0605F_E393F",70,80,2,11,FALSE)
merge2("G0605F_E393F","E577F","G0605F_E393F_E577F",8,15,2,5,FALSE)
merge2("G0605F_E393F_E577F","E783F","G0605F_E393F_E577F_E783F",110,120,4,12,TRUE)
merge2("beewgFor","W158F","beewgFor_W158F",170,180,5,18,TRUE)

#seperate each merged sequence into individual file and build the snp file eventually, for using hisat2 afterward
for (i in 1:4){
  frags <- strsplit(fragName[i],"_")[[1]]
  gene <- genes[i]
  dna <- readDNAStringSet(merge_dir,fragName[i],".fasta"))
  for (j in 1:length(dna)){
    seq <- dna[j]
    name <- names(seq)
    na <- strsplit(name,"_")[[1]]
    sp <- na[1]
    sa <- na[2]
    nbR <- as.integer(na[3:(2+length(frags))])
    ch <- strsplit(as.character(seq),"")[[1]]
    amb <- which(!(ch %in% c("A","T","G","C")))
    if (length(amb)>0){
      if (file.exists(paste0(sampleSpecificRef_dir,gene,"_",name,".snp"))) file.remove(paste0(sampleSpecificRef_dir,gene,"_",name,".snp"))
      for (k in seq_along(amb)){
        al <- IUPAC_CODE_MAP[[ch[amb[k]]]]
        al <- strsplit(al,"")[[1]]
        ch[amb[k]] <- al[1]
        for (a in 2:length(al)){
          cat(k,"\tsingle\t",name,"\t",amb[k]-1,"\t",al[a],"\n",file = paste0(sampleSpecificRef_dir,gene,"_",name,".snp"),append = TRUE)  
        }
      }
      seq <- DNAStringSet(paste0(ch,collapse = ""))
      names(seq) <- name
      writeXStringSet(seq,paste0(sampleSpecificRef_dir,gene,"_",name,".fasta"))
    } else{
      writeXStringSet(seq,paste0(sampleSpecificRef_dir,gene,"_",name,".fasta"))
    }
  }
}
