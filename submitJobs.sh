#step 1 (maximum 1 hour each job)
#The reference fasta file (from sanger sequencing) must be in the folder AGRF_CAGRF15854_B6B68_Reference_fragments/ under the names CO15.fasta etc ...
#The input fastq file must be in the folder AGRF_CAGRF15854_B6B68_Xloose_fragments/ under the names e.g. L.willsi_TA324_M84F.fastq.gz ....
sbatch 1_RemoveContaminants_CountNbReads.sh CO15
sbatch 1_RemoveContaminants_CountNbReads.sh CO13
sbatch 1_RemoveContaminants_CountNbReads.sh EF
sbatch 1_RemoveContaminants_CountNbReads.sh WNT
#step 2 (about 2 hour each job)
sbatch 2_Clustering_RemoveLowPresentSequences.sh M414F
sbatch 2_Clustering_RemoveLowPresentSequences.sh M84F
sbatch 2_Clustering_RemoveLowPresentSequences.sh M202F
sbatch 2_Clustering_RemoveLowPresentSequences.sh M702F
sbatch 2_Clustering_RemoveLowPresentSequences.sh M82F
sbatch 2_Clustering_RemoveLowPresentSequences.sh M414F
sbatch 2_Clustering_RemoveLowPresentSequences.sh G0605F
sbatch 2_Clustering_RemoveLowPresentSequences.sh E393F
sbatch 2_Clustering_RemoveLowPresentSequences.sh E577F
sbatch 2_Clustering_RemoveLowPresentSequences.sh E783F
sbatch 2_Clustering_RemoveLowPresentSequences.sh beewgFor
sbatch 2_Clustering_RemoveLowPresentSequences.sh W158F
#step 3 (<30 minutes)
sbatch runR.sh 3_MergeOverlapSequence.R
#step 4 (maximum 1 hour each job)
sbatch 4_MapToSampleSpecificRef.sh CO15
sbatch 4_MapToSampleSpecificRef.sh CO13
sbatch 4_MapToSampleSpecificRef.sh EF
sbatch 4_MapToSampleSpecificRef.sh WNT
#step5
#The outgroup sequence must be in the folder AGRF_CAGRF15854_B6B68_Reference_fragments/ under the names of CO15_outgroup.fasta etc ...
#Should redefine the path to mafft and RAxML softwares
sbatch 5_MultipleAlignment_Tree.sh CO15
sbatch 5_MultipleAlignment_Tree.sh CO13
sbatch 5_MultipleAlignment_Tree.sh EF
sbatch 5_MultipleAlignment_Tree.sh WNT
#step 6
sbatch 6_MapToSampleSpecificRef_checkReplicates.sh CO13 L.erythrurum  #if you want to check to rebuild consensus sequences of all samples of species L.erythrurum on gene CO13