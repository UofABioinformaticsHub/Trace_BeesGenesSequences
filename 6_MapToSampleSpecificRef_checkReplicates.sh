#!/bin/bash
#SBATCH -p batch                                        # partition (this is the queue your job will be added to)
#SBATCH -N 1                                               # number of nodes
#SBATCH -n 8                                               # number of cores                                          
#SBATCH --time=30:00
#SBATCH --mem=4GB

#This script is a modification of the script in step 4 designed to assemble a common reference sequence between sample/specimen replicates instead of simply assembling one from the major clusters of individual samples.                                                                                             
#The script is only used for samples whose consensus sequences failed to match their replicates yet they were confirmed to be the same using other genes and morphological characters 
#The Gene fragment names used in the manuscript are renamed here as; "CO15_1"= "M414F", "CO15_2" = "M84F", "CO13_1" = "M202F", "CO13_2" = "M702F", "CO13_3" = "M82F", "EF_1" = "G0605F", "EF_2" = "E393F", "EF_3" = "E577F", "EF_4" = "E783F", "WNT_1" = "beewgFor", "WNT_2" = "W158F"

module purge
module load R/3.4.2-foss-2016b
module load HISAT2/2.1.0-foss-2016b
module load seqtk/1.2-foss-2017a
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25
module load BCFtools/1.3.1-GCC-5.3.0-binutils-2.25

INPUT_DIR="AGRF_CAGRF15854_B6B68_MergedXlooseCutadaptDemultiplexed" #Folder contains raw data
ALIGNMENT_DIR="AGRF_CAGRF15854_B6B68_Xloose_alignments_sampleSpecificRef_checkReplicates" #Folder contains alignments to the sample specific reference sequence that is shared between sample replicates
REF_DIR="AGRF_CAGRF15854_B6B68_AGRF_MergedReads_xloose_sampleSpecificRef_checkReplicates" #Folder contains sample specific reference sequences that is shared between sample replicates 
CONSENSUS_DIR="AGRF_CAGRF15854_B6B68_Xloose_Consensus_checkReplicates" #Folder contains consesus sequences generated using the sample specific reference sequence shared between sample replicates. 

mkdir -p ${ALIGNMENT_DIR}
mkdir -p ${CONSENSUS_DIR}
GENE=$1
SPECIES=$2

#select the reference sequences that appears in all replicates of the same species
sed -r "s#gene <-#gene <- \"${GENE}\"#g" checkReplicates.R | sed -r "s#toCheckSpecies <-#toCheckSpecies <- \"${SPECIES}\"#g" > checkReplicates_${GENE}_${SPECIES}.R
R < checkReplicates_${GENE}_${SPECIES}.R --no-save
rm checkReplicates_${GENE}_${SPECIES}.R


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

rm ${CONSENSUS_DIR}/${SPECIES}_consensus.fastq
for FILENAME in $(ls ${REF_DIR}/${GENE}_${SPECIES}_*.fasta); do
    NAME=$(basename $FILENAME)
    NAME=${NAME%.fasta}
    SAMPLE=$(echo ${NAME} | cut -f 3 -d "_")
    
    #build the index for the sample specific reference
    if [ -e ${REF_DIR}/${NAME}.snp ]
    then hisat2-build ${REF_DIR}/${NAME}.fasta ${REF_DIR}/${NAME} --snp ${REF_DIR}/${NAME}.snp
    else hisat2-build ${REF_DIR}/${NAME}.fasta ${REF_DIR}/${NAME}
    fi
    
    #map the raw reads to the new reference for each gene fragment
    FRAG_BAMFILE=""
    for FRAG in ${ALLFRAGS}; do
	    hisat2 -x ${REF_DIR}/${NAME} -U ${INPUT_DIR}/${SAMPLE}_B6B68*${FRAG}.fastq.gz --no-spliced-alignment --ignore-quals --no-softclip  2> ${ALIGNMENT_DIR}/${NAME}_${FRAG}.log | samtools view -Sbh - | samtools sort > ${ALIGNMENT_DIR}/${NAME}_${FRAG}.bam
	    samtools index ${ALIGNMENT_DIR}/${NAME}_${FRAG}.bam
	    FRAG_BAMFILE=$(echo ${FRAG_BAMFILE} "${ALIGNMENT_DIR}/${NAME}_${FRAG}.bam")
    done
    
    #merge the alignments of all gene fragments into one bam file
    samtools merge -f ${ALIGNMENT_DIR}/${NAME}.bam ${FRAG_BAMFILE}
    samtools index ${ALIGNMENT_DIR}/${NAME}.bam
    
    #build the consensus sequnces from the alignments
    samtools mpileup -B -d 100000 -u ${ALIGNMENT_DIR}/${NAME}.bam | bcftools call -c -M vcfutils.pl vcf2fq - > ${CONSENSUS_DIR}/${NAME}_consensus.fastq
    nb=$(wc -l ${CONSENSUS_DIR}/${NAME}_consensus.fastq | cut -f1 -d " ")
    if [ $nb -ne 2 ] 
    then 
	cat ${CONSENSUS_DIR}/${NAME}_consensus.fastq >> ${CONSENSUS_DIR}/${GENE}_${SPECIES}_consensus.fastq #append consensus sequence of all samples into one file
    fi
done

seqtk seq -a ${CONSENSUS_DIR}/${GENE}_${SPECIES}_consensus.fastq > ${CONSENSUS_DIR}/${GENE}_${SPECIES}_consensus.fasta #convert fastq to fasta file

