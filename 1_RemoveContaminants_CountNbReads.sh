#!/bin/bash
#SBATCH -p batch                                        # partition (this is the queue your job will be added to)
#SBATCH -N 1                                               # number of nodes
#SBATCH -n 8                                               # number of cores                                          
#SBATCH --time=2:00:00
#SBATCH --mem=4GB

#This script maps raw reads to reference sequences and retains only those reads that are 84% similar to the reference sequence  
module purge
module load HISAT2/2.1.0-foss-2016b
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25

#Set up
GENE=$1  ## ex: CO13
REF_DIR=AGRF_CAGRF15854_B6B68_Reference_fragments  ## reference folder
INPUT_DIR=AGRF_CAGRF15854_B6B68_MergedXlooseCutadaptDemultiplexed   ## input folder contains fastq files
ALIGNMENT_DIR=AGRF_CAGRF15854_B6B68_Merged_Genefragment_trimmed_alignments_fragments  ## folder contains alignments
COUNT_DIR=AGRF_CAGRF15854_B6B68_MergedReads_xloose_readCounts ##Folder contains read counts for each sample and gene fragment
mkdir -p ${ALIGNMENT_DIR}

#Build indexes for the reference sequences   
if [ -e ${REF_DIR}/${GENE}.snp ]
then hisat2-build ${REF_DIR}/${GENE}.fasta ${REF_DIR}/${GENE} --snp ${REF_DIR}/${GENE}.snp
else hisat2-build ${REF_DIR}/${GENE}.fasta ${REF_DIR}/${GENE}
fi

if [ ${GENE} == "CO13" ]
then
    FRAG="M202F M702F M82F"
fi
if [ ${GENE} == "CO15" ]
then
    FRAG="M414F M84F"
fi
if [ ${GENE} == "EF" ]
then
    FRAG="G0605F E393F E577F E783F"
fi
if [ ${GENE} == "WNT" ]
then
    FRAG="beewgFor W158F"
fi

for f in ${FRAG}; do
  rm -f ${COUNT_DIR}/${f}.txt
  rm -f ${COUNT_DIR}/${f}.nb
  for FILE in ${INPUT_DIR}/*_B6B68_*-*-${f}.fastq.gz; do
  #Map the raw reads of each fragment to the reference sequences and count number of reads passing filter
      R=$(basename $FILE)     
      R=${R%.fastq.gz}
      hisat2 -x ${REF_DIR}/${GENE} -U ${INPUT_DIR}/${R}.fastq.gz --sp 6,2 --score-min L,0,-1  --no-spliced-alignment --ignore-quals  2> ${ALIGNMENT_DIR}/${R}.log | samtools view -Sbh - > ${ALIGNMENT_DIR}/${R}.bam
      samtools view -F 260 ${ALIGNMENT_DIR}/${R}.bam | cut -f 10 | awk '{if ((length($1) > 250)) print;}'| sort | uniq -c | sort -nr | sed -r 's#^ +##g' | sed -r 's#\t# #g' > ${COUNT_DIR}/${R}.txt 
      echo ${R} ${f}
      cat ${COUNT_DIR}/${R}.txt >> ${COUNT_DIR}/${f}.txt
      wc -l ${COUNT_DIR}/${R}.txt >> ${COUNT_DIR}/${f}.nb
   done
done
