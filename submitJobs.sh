#The scripts can be submitted to a high performance computing facility in the order shown below.
sbatch 0_Script_For_Analysing_Rawdata.sh
#step 1 (maximum 1 hour each job)
#The reference fasta file (from sanger sequencing) must be in the folder AGRF_CAGRF15854_B6B68_Reference_fragments/ under the names CO15.fasta etc ...
#The input fastq file must be in the folder AGRF_CAGRF15854_B6B68_MergedXlooseCutadaptDemultiplexed under the names e.g. TA020.1_B6B68_GGACTCCT-CTCTCTAT-E393F.fastq.gz.. etc
sbatch 1_RemoveContaminants_CountNbReads.sh CO15
sbatch 1_RemoveContaminants_CountNbReads.sh CO13
sbatch 1_RemoveContaminants_CountNbReads.sh EF
sbatch 1_RemoveContaminants_CountNbReads.sh WNT
#step 2 (about 2 hour each job)
sbatch 2_Cluster_RemoveSmallClusters.sh M414F
sbatch 2_Cluster_RemoveSmallClusters.sh M84F
sbatch 2_Cluster_RemoveSmallClusters.sh M202F
sbatch 2_Cluster_RemoveSmallClusters.sh M702F
sbatch 2_Cluster_RemoveSmallClusters.sh M82F
sbatch 2_Cluster_RemoveSmallClusters.sh M414F
sbatch 2_Cluster_RemoveSmallClusters.sh G0605F
sbatch 2_Cluster_RemoveSmallClusters.sh E393F
sbatch 2_Cluster_RemoveSmallClusters.sh E577F
sbatch 2_Cluster_RemoveSmallClusters.sh E783F
sbatch 2_Cluster_RemoveSmallClusters.sh beewgFor
sbatch 2_Cluster_RemoveSmallClusters.sh W158F
#step 3 (<30 minutes)
sbatch runR.sh 3_MergeOverlapSequence.R
#step 4 (maximum 1 hour each job)
sbatch 4_MapToSampleSpecificRef.sh CO15
sbatch 4_MapToSampleSpecificRef.sh CO13
sbatch 4_MapToSampleSpecificRef.sh EF
sbatch 4_MapToSampleSpecificRef.sh WNT
#step5
#The outgroup sequence must be in the folder AGRF_CAGRF15854_B6B68_Reference_fragments/ under the names of CO15_outgroup.fasta etc ...
sbatch 5_MultipleAlignment_Tree.sh CO15
sbatch 5_MultipleAlignment_Tree.sh CO13
sbatch 5_MultipleAlignment_Tree.sh EF
sbatch 5_MultipleAlignment_Tree.sh WNT
#step 6
#Used to rescue samples whose consensus sequences do not match between replicates if the identification of the replicates is confirmed to be accurate
sbatch 6_MapToSampleSpecificRef_checkReplicates.sh CO13 L.erythrurum  