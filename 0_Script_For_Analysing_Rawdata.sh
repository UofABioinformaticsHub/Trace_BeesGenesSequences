#!/bin/bash
#SBATCH -p batch                                        # partition (this is the queue your job will be added to)
#SBATCH -N 1                                               # number of nodes
#SBATCH -n 8                                              # number of cores
#SBATCH --time=20:00:00
#SBATCH --mem=8GB

#Replace paths to input files appropriately 
#Uses BBMERGE to combine R1 and R2 of the same samples and outputs compressed files

module load BBMap/35.92-GCC-5.3.0-binutils-2.25-Java-1.8.0_71 

NAMES=$(ls AGRF_CAGRF15854_B6B68_RawData/*B6B68*L001_R1.fastq.gz | sed 's#1.fastq.gz##g')
for NAME in ${NAMES}; do
    bbmerge.sh -xloose in1=${NAME}1.fastq.gz in2=${NAME}2.fastq.gz out=./AGRF_CAGRF15854_B6B68_MergedReads_xloose/${NAME}.merged.fastq.gz outu1=./AGRF_CAGRF15854_B6B68unmerged1_xloose/${NAME}1.unmerged.fastq.gz outu2=./AGRF_CAGRF15854_B6B68unmerged2_xloose/${NAME}2.unmerged.fastq.gz 
done

#Loops through BBmerged files and trims the forward anchored primers while demultiplexing them into individual files by sequenced fragment..
#A fasta file with the Forward primers to be trimmed should be provided in the format described on the cutadapt website 

module load cutadapt/1.9.1-foss-2016b-Python-2.7.13

for f in $(ls AGRF_CAGRF15854_B6B68_MergedReads_xloose/*.merged.fastq.gz); do
    base=$(basename $f)
    base=${base%_L001_R.merged.fastq.gz};
    outfile=AGRF_CAGRF15854_B6B68_MergedXlooseCutadapt/${base};
    untrimfile=AGRF_CAGRF15854_B6B68_MergedXlooseCutadapt/untrimmed_${base};
    cutadapt -g file:AGRF_CAGRF15854_B6B68_F.fasta -o ${outfile}-{name}.fastq.gz --untrimmed-output ${untrimfile}.fastq.gz $f;
done

#Loops through demultiplexed files of each gene fragment (labeled by primer name in list 'GENEFRAG' trimming anchored specific reverse degenerate primers (RPRIMER)

GENEFRAG=(M414F M84F M202F M702F M82F beewgFor W158F G0605F E393F E577F E783F)
RPRIMER=(GAACARTWTAYCCHCCHYTATC TGATTTTTTGGHCAYCCWGAAGTWTA GGDGGNHTWACWGGWATYAT CCWCGWCGWTAYTCWGAYTAYCC TAATATGGCAGATTAGTCGATTGGA CACGGYAGACAGTGYAACGA TGYACATTCCAYTGGTGYTGCGNAGT CTCCGAAGCYCGATTYGAAG GGCTCTCCGTCTTCCYCTTCAGG CCCGTTGGTCGTGTCGAAACT CCACCTAAAGGTGCTGCTGATT)
for i in {1..11}; do
    for f in $(ls AGRF_CAGRF15854_B6B68_MergedXlooseCutadapt/*${GENEFRAG[$i]}.fastq.gz); do
	base=$(basename $f)
	base=${base%.fastq.gz};
	outfile=AGRF_CAGRF15854_B6B68_MergedXlooseCutadaptDemultiplexed/${base};
	untrimfile=AGRF_CAGRF15854_B6B68_MergedXlooseCutadaptDemultiplexed/untrimmed_${base};
	cutadapt -a ${RPRIMER[$i]}$  -o ${outfile}.fastq.gz --untrimmed-output ${untrimfile}.fastq.gz $f;
    done
done

