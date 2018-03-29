source("utils.R")
gene <- 
toCheckSpecies <- 
  
if (gene=="CO15"){
  frags <- c("M414F","M84F")
}
if (gene=="CO13"){
  frags <- c("M202F","M702F","M82F")
}
if (gene=="EF"){
  frags <- c("G0605F","E393F","E577F","E783F")
}
if (gene=="WNT"){
  frags <- c("beewgFor","W158F")
}

for (frag in frags){
    sequences <- readDNAStringSet(paste0(merge_dir,frag,".fasta"))
    names <- names(sequences)
    ss <- strsplit(names,":")
    ss <- sapply(ss, function(s){s[length(s)]})
    ss <- strsplit(ss,"_")
    species <- sapply(ss,function(i){i[2]})
    writeXStringSet(sequences[species==toCheckSpecies],paste0(merge_dir,frag,"_",toCheckSpecies,".fasta"))
}

if (gene=="CO15"){
  merge3(toCheckSpecies,"M414F","M84F","M414F_M84F",7,20,2,3,FALSE)
}
if (gene=="CO13"){
  merge3(toCheckSpecies,"M202F","M702F","M202F_M702F",60,70,2,10,FALSE)
  merge3(toCheckSpecies,"M202F_M702F","M82F","M202F_M702F_M82F",20,30,2,5,TRUE)
}
if (gene=="EF"){
  merge3(toCheckSpecies,"G0605F","E393F","G0605F_E393F",70,80,2,11,FALSE)
  merge3(toCheckSpecies,"G0605F_E393F","E577F","G0605F_E393F_E577F",8,15,2,5,FALSE)
  merge3(toCheckSpecies,"G0605F_E393F_E577F","E783F","G0605F_E393F_E577F_E783F",110,120,4,12,TRUE)
}
if (gene=="WNT"){
  merge3(toCheckSpecies,"beewgFor","W158F","beewgFor_W158F",170,180,5,18,TRUE)
}

#seperate each merged sequence into individual file and build the snp file eventually, for using hisat2 afterward
fragName <- paste0(frags,collapse = "_")
dna <- readDNAStringSet(paste0(merge_dir,fragName,"_",toCheckSpecies,".fasta"))
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
    if (file.exists(paste0(sampleSpecificRef_checkReplicates_dir,gene,"_",name,".snp"))) file.remove(paste0(sampleSpecificRef_checkReplicates_dir,gene,"_",name,".snp"))
    for (k in seq_along(amb)){
      al <- IUPAC_CODE_MAP[[ch[amb[k]]]]
      al <- strsplit(al,"")[[1]]
      ch[amb[k]] <- al[1]
      for (a in 2:length(al)){
        cat(k,"\tsingle\t",name,"\t",amb[k]-1,"\t",al[a],"\n",file = paste0(sampleSpecificRef_checkReplicates_dir,gene,"_",name,".snp"),append = TRUE)  
      }
    }
    seq <- DNAStringSet(paste0(ch,collapse = ""))
    names(seq) <- name
    writeXStringSet(seq,paste0(sampleSpecificRef_checkReplicates_dir,gene,"_",name,".fasta"))
  } else{
    writeXStringSet(seq,paste0(sampleSpecificRef_checkReplicates_dir,gene,"_",name,".fasta"))
  }
}
