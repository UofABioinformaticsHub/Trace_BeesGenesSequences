#!/bin/bash
#SBATCH -p batch                                        # partition (this is the queue your job will be added to)
#SBATCH -N 1                                               # number of nodes
#SBATCH -n 8                                               # number of cores                                          
#SBATCH --time=2:00:00
#SBATCH --mem=4GB


#This script maps raw reads to their sample specific reference sequnces and generates consesus sequences for each sample and gene
#The Gene fragment names used in the manuscript are renamed here as; "CO15_1"= "M414F", "CO15_2" = "M84F", "CO13_1" = "M202F", "CO13_2" = "M702F", "CO13_3" = "M82F", "EF_1" = "G0605F", "EF_2" = "E393F", "EF_3" = "E577F", "EF_4" = "E783F", "WNT_1" = "beewgFor", "WNT_2" = "W158F"

module purge
module load R/3.4.4-foss-2016b
module load HISAT2/2.1.0-foss-2016b
module load seqtk/1.2-foss-2017a
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25
module load BCFtools/1.3.1-GCC-5.3.0-binutils-2.25


input_dir="AGRF_CAGRF15854_B6B68_MergedXlooseCutadaptDemultiplexed"## input folder contains raw data fastq files  
alignment_dir="AGRF_CAGRF15854_B6B68_Merged_Genefragment_trimmed_alignments_sampleSpecificRef"## Folder contains alignments to sample specifc reference sequences
ref_dir="AGRF_CAGRF15854_B6B68_MergedReads_xloose_sampleSpecificRef"##Folder contains reference sequences for individual samples
CONSENSUS_DIR="AGRF_CAGRF15854_B6B68_Merged_Genefragment_trimmed_consensus_sampleSpecificRef"# Folder contains assembled consesus sequences for each sample
mkdir -p ${alignment_dir}
mkdir -p ${CONSENSUS_DIR}
gene=$1
if [ ${gene} == "CO13" ]
then
    FRAG="M202F M702F M82F"
fi
if [ ${gene} == "CO15" ]
then
    FRAG="M414F M84F"
fi
if [ ${gene} == "EF" ]
then
    FRAG="G0605F E393F E577F E783F"
fi
if [ ${gene} == "WNT" ]
then
    FRAG="beewgFor W158F"
fi


rm ${CONSENSUS_DIR}/${gene}_consensus.fastq
for name in $(ls ${ref_dir}/${gene}_*TA*.fasta); do
    n=$(basename $name)
    n=${n%.fasta}
    sa=$(echo $n | cut -f 2 -d "_")

    #build the hisat2 index for the reference sequence
    if [ -e ${ref_dir}/${n}.snp ]
    then hisat2-build ${ref_dir}/${n}.fasta ${ref_dir}/${n} --snp ${ref_dir}/${n}.snp
    else hisat2-build ${ref_dir}/${n}.fasta ${ref_dir}/${n}
    fi
    
    #map the reads of each gene fragment to the reference
    FRAG_BAMFILE=""
    for f in ${FRAG}; do
	hisat2 -x ${ref_dir}/${n} -U ${input_dir}/${sa}_B6B68_*-*-${f}.fastq.gz --no-spliced-alignment --ignore-quals --no-softclip  2> ${alignment_dir}/${n}_${f}.log | samtools view -Sbh - | samtools sort > ${alignment_dir}/${n}_${f}.bam
	samtools index ${alignment_dir}/${n}_${f}.bam
	FRAG_BAMFILE=$(echo ${FRAG_BAMFILE} "${alignment_dir}/${n}_${f}.bam")
    done
    #merge the alignments of all gene fragments into one file
    samtools merge -f ${alignment_dir}/${n}.bam ${FRAG_BAMFILE}
    samtools index ${alignment_dir}/${n}.bam

    #generate the consensus sequence from the merged alignments
    samtools mpileup -B -d 100000 -u ${alignment_dir}/${n}.bam | bcftools call -c -M | vcfutils.pl vcf2fq - > ${CONSENSUS_DIR}/${n}_consensus.fastq
    nb=$(wc -l ${CONSENSUS_DIR}/${n}_consensus.fastq | cut -f1 -d " ")
    if [ $nb -ne 2 ] 
    then 
	cat ${CONSENSUS_DIR}/${n}_consensus.fastq >> ${CONSENSUS_DIR}/${gene}_consensus.fastq
    fi
done

#convert the fastq consensus sequence file into fasta file
seqtk seq -a ${CONSENSUS_DIR}/${gene}_consensus.fastq > ${CONSENSUS_DIR}/${gene}_consensus.fasta 

# Rename: add species names into consensus sequence names
# The species name is defined in the file 'Species.txt' 

sed -r "s#gene <-#gene <- \"${gene}\"#g" addSpeciesToConsensusName.R > addSpeciesToConsensusName_${gene}.R
R < addSpeciesToConsensusName_${gene}.R --no-save
rm addSpeciesToConsensusName_${gene}.R
