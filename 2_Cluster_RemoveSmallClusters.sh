#!/bin/bash
#SBATCH -p batch                                        # partition (this is the queue your job will be added to) 
#SBATCH -N 1                                               # number of nodes 
#SBATCH -n 1                                              # number of cores 
#SBATCH --time=2:30:00                                    # time allocation
#SBATCH --mem=10GB 

#This script uses an R script to cluster all reads at 98% and removes any clusters having a low number of reads (less that 1% of the total number of reads in a sample and less than 3 reads in total). 
module load R/3.4.2-foss-2016b

FRAG=$1
sed -r "s#frag <-#frag <- \"${FRAG}\"#g" Cluster_RemoveSmallClusters.R > Cluster_RemoveSmallClusters_${FRAG}.R
R < Cluster_RemoveSmallClusters_${FRAG}.R --no-save
rm Cluster_RemoveSmallClusters_${FRAG}.R
