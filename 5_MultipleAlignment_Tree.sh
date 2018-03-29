#!/bin/bash
#SBATCH -p batch                                        # partition (this is the queue your job will be added to)
#SBATCH -N 1                                               # number of nodes
#SBATCH -n 8                                               # number of cores                                          
#SBATCH --time=2:00:00
#SBATCH --mem=4GB

GENE=$1
CONSENSUS_DIR="AGRF_CAGRF15854_B6B68_Xloose_Consensus"
REF_DIR="AGRF_CAGRF15854_B6B68_Reference_fragments"
TREE_DIR="$PWD/AGRF_CAGRF15854_B6B68_Xloose_RAxML"
mkdir -p ${TREE_DIR}

MAFFT=/fast/users/a1221455/mafft/bin/mafft
RAXML=/fast/users/a1221455/standard-RAxML/raxmlHPC

#add the outgroup sequence into the consensus sequences file. The outgroup sequence must be in the foler ${REF_DIR}
cat ${CONSENSUS_DIR}/${GENE}_consensus.fasta ${REF_DIR}/${GENE}_outgroup.fasta > ${CONSENSUS_DIR}/${GENE}_consensus_outgroup.fasta 

#multiple alignment
${MAFFT} ${CONSENSUS_DIR}/${GENE}_consensus_outgroup.fasta > ${TREE_DIR}/${GENE}_consensus_outgroup.aln

#building trees using RAXML
rm ${TREE_DIR}/RAxML*${GENE}
${RAXML} -p 12345 -s ${TREE_DIR}/${GENE}_consensus_outgroup.aln -m GTRGAMMAI -o H.megastigmus -n ${GENE} -w ${TREE_DIR}

