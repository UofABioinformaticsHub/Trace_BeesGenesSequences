#!/bin/bash
#SBATCH -p batch                                        # partition (this is the queue your job will be added to) 
#SBATCH -N 1                                               # number of nodes 
#SBATCH -n 1                                              # number of cores 
#SBATCH --time=2:30:00                                    # time allocation
#SBATCH --mem=10GB 
module load R/3.4.0-foss-2016b

FRAG=$1
sed -r "s#frag <-#frag <- \"${FRAG}\"#g" Clustering_RemoveLowPresentSequences.R > Clustering_RemoveLowPresentSequences_${FRAG}.R
R < Clustering_RemoveLowPresentSequences_${FRAG}.R --no-save
rm Clustering_RemoveLowPresentSequences_${FRAG}.R
