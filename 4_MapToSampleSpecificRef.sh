#!/bin/bash
#SBATCH -p batch                                        # partition (this is the queue your job will be added to)
#SBATCH -N 1                                               # number of nodes
#SBATCH -n 8                                               # number of cores                                          
#SBATCH --time=2:00:00
#SBATCH --mem=4GB

module purge
module load HISAT2/2.1.0-foss-2016b
module load seqtk/1.2-foss-2017a
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25
module load BCFtools/1.3.1-GCC-5.3.0-binutils-2.25

INPUT_FIR="AGRF_CAGRF15854_B6B68_Xloose_fragments"
REF_DIR="AGRF_CAGRF15854_B6B68_Xloose_sampleSpecificRef"
ALIGNMENT_DIR="AGRF_CAGRF15854_B6B68_Xloose_alignments_sampleSpecificRef"
CONSENSUS_DIR="AGRF_CAGRF15854_B6B68_Xloose_Consensus"

mkdir -p ${ALIGNMENT_DIR}
mkdir -p ${CONSENSUS_DIR}
GENE=$1
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

rm ${CONSENSUS_DIR}/${GENE}_consensus.fastq
for FILENAME in $(ls ${REF_DIR}/${GENE}_*TA*.fasta); do
    NAME=$(basename $FILENAME)
    NAME=${NAME%.fasta}
    SPECIES=$(echo $NAME | cut -f 2 -d "_")
    SAMPLE=$(echo $NAME | cut -f 3 -d "_")
    
    #build index for the sample specific reference
    if [ -e ${REF_DIR}/${NAME}.snp ]
      then hisat2-build ${REF_DIR}/${NAME}.fasta ${REF_DIR}/${NAME} --snp ${REF_DIR}/${NAME}.snp
      else hisat2-build ${REF_DIR}/${NAME}.fasta ${REF_DIR}/${NAME}
    fi
    FRAG_BAMFILE=""
    
    #map the raw reads to the sample specific reference for each gene fragment
    for FRAG in ${ALLFRAGS}; do
	    hisat2 -x ${REF_DIR}/${NAME} -U ${INPUT_FIR}/${SPECIES}_${SAMPLE}_${FRAG}.fastq.gz --no-spliced-alignment --ignore-quals --no-softclip  2> ${ALIGNMENT_DIR}/${NAME}_${FRAG}.log | samtools view -Sbh - | samtools sort > ${ALIGNMENT_DIR}/${NAME}_${FRAG}.bam
	    samtools index ${ALIGNMENT_DIR}/${NAME}_${FRAG}.bam
	    FRAG_BAMFILE=$(echo ${FRAG_BAMFILE} "${ALIGNMENT_DIR}/${NAME}_${FRAG}.bam")
    done
    #merge the alignments of all gene fragment into one file
    samtools merge -f ${ALIGNMENT_DIR}/${NAME}.bam ${FRAG_BAMFILE}
    samtools index ${ALIGNMENT_DIR}/${NAME}.bam
    
    #build the consensus sequence from the alignment
    samtools mpileup -B -d 100000 -u ${ALIGNMENT_DIR}/${NAME}.bam | bcftools call -c -M > ${CONSENSUS_DIR}/${NAME}.bcf
    samtools mpileup -B -d 100000  ${ALIGNMENT_DIR}/${NAME}.bam > ${CONSENSUS_DIR}/${NAME}.coverage 
    vcfutils.pl vcf2fq ${CONSENSUS_DIR}/${NAME}.bcf > ${CONSENSUS_DIR}/${NAME}_consensus.fastq
    nb=$(wc -l ${CONSENSUS_DIR}/${NAME}_consensus.fastq | cut -f1 -d " ")
    if [ $nb -ne 2 ] 
    then 
	     cat ${CONSENSUS_DIR}/${NAME}_consensus.fastq >> ${CONSENSUS_DIR}/${GENE}_consensus.fastq #append consensus sequence of all samples into one file
    fi
done

seqtk seq -a ${CONSENSUS_DIR}/${GENE}_consensus.fastq > ${CONSENSUS_DIR}/${GENE}_consensus.fasta #convert fastq to fasta file


