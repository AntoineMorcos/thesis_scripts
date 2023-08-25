#!/bin/bash -l

#SBATCH --output=slurm-hisat2-%u-%j.out  ## output file name
#SBATCH --job-name=hisat2 ## name of job we'll see in the queue
#SBATCH --partition=brc
#SBATCH --nodes=1 # use 1 node (unless the code is multi-node parallelized)
#SBATCH --mem=8G ## 8GB of RAM


echo "start: " `date`

# load the module
module load apps/hisat2/2.1.0-python3.7.3

# Set the path to the genome index files indicating only the prefix
genome_index="/scratch/groups/rosetree/RNAseq_2020/hg38_tran/genome_tran"
# in the directory Hisat2_Index_GRCh38 there are 8 files: genome_tran.1.ht2 to genome_tran.8.ht2

# Set the input read file paths - paired end reads
fwd_file="/scratch/groups/rosetree/RNAseq_2020/raw_data/all_files/L1_ID1_1.fq.gz"
rev_file="/scratch/groups/rosetree/RNAseq_2020/raw_data/all_files/L1_ID1_2.fq.gz"

# Set the output file path
alignment_output_file="/scratch/groups/rosetree/RNAseq_2020/hisat2_alignment/L1_ID1_hisat2.sam"

# Summary file path
summary_file="/scratch/groups/rosetree/RNAseq_2020/hisat2_alignment/summary_files/L1_ID1_hisat2_summary.txt"



# Run HISAT2
hisat2 -p 12 -x $genome_index -1 $fwd_file -2 $rev_file -S $alignment_output_file --summary-file $summary_file
# Usage: hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]
# -p is number of alignment threads to launch
# -x is the Index filename prefix (minus trailing .X.ht2)
# -1 Files with #1 mates, paired with files in <m2>. Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
# -2 Files with #2 mates, paired with files in <m1>. Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
# -S File for SAM output (default: stdout)

echo "end: " `date`