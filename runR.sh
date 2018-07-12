#!/bin/bash
#SBATCH -p batch                                        # partition (this is the queue your job will be added to) 
#SBATCH -N 1                                               # number of nodes 
#SBATCH -n 1                                              # number of cores 
#SBATCH --time=3:00:00                                    # time allocation
#SBATCH --mem=10GB 
module load R/3.4.2-foss-2016b
R < $1 --no-save
