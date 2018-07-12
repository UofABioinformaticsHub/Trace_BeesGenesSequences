#!/bin/bash
#SBATCH -p batch                                        # partition (this is the queue your job will be added to)
#SBATCH -N 1                                               # number of nodes
#SBATCH -n 8                                               # number of cores                                          
#SBATCH --time=2:00:00
#SBATCH --mem=4GB

#This script creates RAxML phylogenetic trees using the consensus sequences from step 4

module load MUSCLE/3.8.31

CONSENSUS_DIR="AGRF_CAGRF15854_B6B68_Merged_Genefragment_trimmed_consensus_sampleSpecificRef" #Folder contains consesus sequences for each sample and gene
ref_dir="AGRF_CAGRF15854_B6B68_Reference_fragments" #Folder contains sequences of the outgroup  
TREE_DIR="${PWD}/AGRF_CAGRF15854_B6B68_xloose_RAxML" #Folder contains RAXML trees generated 
gene=$1
mkdir -p ${TREE_DIR}
cat ${CONSENSUS_DIR}/${gene}_consensus_speciesName.fasta ${ref_dir}/${gene}_outgroup.fasta > ${CONSENSUS_DIR}/${gene}_consensus_outgroup.fasta 
muscle -in ${CONSENSUS_DIR}/${gene}_consensus_outgroup.fasta -out ${CONSENSUS_DIR}/${gene}_consensus_outgroup.aln
rm ${TREE_DIR}/RAxML*${gene}
/fast/users/a1648633/standard-RAxML/raxmlHPC -p 12345 -s ${CONSENSUS_DIR}/${gene}_consensus_outgroup.aln -m GTRGAMMAI -o H.megastigmus -n ${gene} -w ${TREE_DIR}

