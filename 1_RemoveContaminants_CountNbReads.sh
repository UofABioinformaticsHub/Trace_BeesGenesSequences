#!/bin/bash
#SBATCH -p batch                                        # partition (this is the queue your job will be added to)
#SBATCH -N 1                                               # number of nodes
#SBATCH -n 8                                               # number of cores                                          
#SBATCH --time=1:00:00
#SBATCH --mem=4GB

module purge
module load HISAT2/2.1.0-foss-2016b
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25

#SET UP
GENE=$1  ## ex: CO13
REF_DIR=AGRF_CAGRF15854_B6B68_Reference_fragments  ## reference folder, contains the reference of each gene (from Sanger sequencing)
INPUT_DIR=AGRF_CAGRF15854_B6B68_Xloose_fragments   ## input folder contains fastq files e.g. L.willsi_TA324_M84F.fastq.gz
ALIGNMENT_DIR=AGRF_CAGRF15854_B6B68_Xloose_alignments  ## folder contains alignments bam file
COUNT_DIR=AGRF_CAGRF15854_B6B68_Xloose_readCounts ## folder contains the read count files
mkdir -p ${ALIGNMENT_DIR}
mkdir -p ${COUNT_DIR}

#BUILD THE INDEX OF REFERENCE   
if [ -e ${REF_DIR}/${GENE}.snp ]
then hisat2-build ${REF_DIR}/${GENE}.fasta ${REF_DIR}/${GENE} --snp ${REF_DIR}/${GENE}.snp
else hisat2-build ${REF_DIR}/${GENE}.fasta ${REF_DIR}/${GENE}
fi

if [ ${GENE} == "CO13" ]
then
    ALLFRAGS="M202F M702F M82F"
fi
if [ ${GENE} == "CO15" ]
then
    ALLFRAGS="M414F M84F"
fi
if [ ${GENE} == "EF" ]
then
    ALLFRAGS="G0605F E393F E577F E783F"
fi
if [ ${GENE} == "WNT" ]
then
    ALLFRAGS="beewgFor W158F"
fi

#Define all available species
ALLSPECIES=$((ls ${INPUT_DIR}/*_*_E577F.fastq.gz) | sed -r "s#${INPUT_DIR}\/##g" | sed -r "s#_.+_E577F.fastq.gz##g" | sort | uniq)

for FRAG in ${ALLFRAGS}; do
  rm -f ${COUNT_DIR}/${FRAG}.txt
  rm -f ${COUNT_DIR}/${FRAG}.nb
  for SPECIES in ${ALLSPECIES}; do
      SAMPLES=$(ls ${INPUT_DIR}/${SPECIES}_*_${FRAG}.fastq.gz)
      rm  -f ${COUNT_DIR}/${SPECIES}_${FRAG}.txt
      rm -f ${COUNT_DIR}/${SPECIES}_${FRAG}.nb
      for SAMPLE in ${SAMPLES}; do
      #MAP THE RAW READS OF EACH FRAGMENT TO REFERENCE
	      SAMPLE=$(basename $SAMPLE)     
	      SAMPLE=${SAMPLE%.fastq.gz}
	      hisat2 -x ${REF_DIR}/${GENE} -U ${INPUT_DIR}/${SAMPLE}.fastq.gz --sp 6,2 --score-min L,0,-1  --no-spliced-alignment --ignore-quals  2> ${ALIGNMENT_DIR}/${SAMPLE}.log | samtools view -Sbh - > ${ALIGNMENT_DIR}/${SAMPLE}.bam
	      samtools view -F 260 ${ALIGNMENT_DIR}/${SAMPLE}.bam | cut -f 10 | awk '{if ((length($1) > 250)) print;}'| sort | uniq -c | sort -nr | sed -r 's#^ +##g' | sed -r 's#\t# #g' > ${COUNT_DIR}/${SAMPLE}.txt 
	      cat ${COUNT_DIR}/${SAMPLE}.txt >> ${COUNT_DIR}/${SPECIES}_${FRAG}.txt
	      wc -l ${COUNT_DIR}/${SAMPLE}.txt >> ${COUNT_DIR}/${SPECIES}_${FRAG}.nb
      done
      cat ${COUNT_DIR}/${SPECIES}_${FRAG}.txt >> ${COUNT_DIR}/${FRAG}.txt
          wc -l ${COUNT_DIR}/${SPECIES}_${FRAG}.txt >> ${COUNT_DIR}/${FRAG}.nb
    done
done
