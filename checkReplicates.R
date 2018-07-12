source("utils.R")
gene <- #to put the gene name to anlysis
toCheckSpecies <- #to put the species name to check 
  
#Define the gene fragment names for each gene
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

#Read the species names given for each sample
mapSpecies <- as.data.frame(read.table("Species.txt"))
colnames(mapSpecies) <- c("Species","Sample")

#Extract the sequences of samples to check
for (frag in frags){
  sequences <- readDNAStringSet(paste0(merge_dir,frag,".fasta"))
  names <- names(sequences)
  samples <- strsplit(names,":") %>% 
    sapply(function(s){s[length(s)]}) %>% 
    strsplit("_") %>% 
    sapply(function(i){i[2]})
  idSamples <- data.frame(Id = 1:length(samples), Sample = samples)
  idSamples <- left_join(idSamples,mapSpecies,by="Sample")
  toCheckSamples <- which(idSamples$Species==toCheckSpecies)
  writeXStringSet(sequences[toCheckSamples],paste0(merge_dir,frag,"_",toCheckSpecies,".fasta"))
}

#Merge the overlapped sequences of fragments coming from the same sample and select the most
#consistent one for each sample
if (gene=="CO15"){
  merge3(toCheckSpecies,"M414F","M84F","M414F_M84F",7,20,2,3,TRUE)
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

#Seperate each merged sequence into individual files and build the snp file for using hisat2
fragName <- paste0(frags,collapse = "_")
dna <- readDNAStringSet(paste0(merge_dir,fragName,"_",toCheckSpecies,".fasta"))
for (j in 1:length(dna)){
  seq <- dna[j]
  name <- names(seq)
  na <- strsplit(name,"_")[[1]]
  sa <- na[1]
  nbR <- as.integer(na[2:(1+length(frags))])
  ch <- strsplit(as.character(seq),"")[[1]]
  amb <- which(!(ch %in% c("A","T","G","C")))
  if (length(amb)>0){
    if (file.exists(paste0(sampleSpecificRef_checkReplicates_dir,gene,"_",toCheckSpecies,"_",name,".snp"))){ 
      file.remove(paste0(sampleSpecificRef_checkReplicates_dir,gene,"_",toCheckSpecies,"_",name,".snp"))
    }
    for (k in seq_along(amb)){
      al <- IUPAC_CODE_MAP[[ch[amb[k]]]]
      al <- strsplit(al,"")[[1]]
      ch[amb[k]] <- al[1]
      for (a in 2:length(al)){
        cat(k,"\tsingle\t",name,"\t",amb[k]-1,"\t",al[a],"\n",
            file = paste0(sampleSpecificRef_checkReplicates_dir,gene,"_",toCheckSpecies,"_",name,".snp"),append = TRUE)  
      }
    }
    seq <- DNAStringSet(paste0(ch,collapse = ""))
    names(seq) <- name
    writeXStringSet(seq,paste0(sampleSpecificRef_checkReplicates_dir,gene,"_",toCheckSpecies,"_",name,".fasta"))
  } else{
    writeXStringSet(seq,paste0(sampleSpecificRef_checkReplicates_dir,gene,"_",toCheckSpecies,"_",name,".fasta"))
  }
}