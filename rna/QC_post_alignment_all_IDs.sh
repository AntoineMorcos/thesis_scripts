#!/bin/bash -l

#SBATCH --output=slurm-QC-loop-after-alignment-%u-%j.out  ## output file name
#SBATCH --job-name=QC ## name of job we'll see in the queue
#SBATCH --partition=brc
#SBATCH --nodes=1 # use 1 node (unless the code is multi-node parallelized)
#SBATCH --mem=8G ## 8GB of RAM



########################################################################
#                         Basic set up                                 #
########################################################################

# print start time to slurm output file
echo "start: " `date`

# load the modules
module load apps/samtools/1.9.0-singularity
module load apps/fastqc/0.11.8



########################################################################
#                     1. Conversion of SAM to BAM                      #
########################################################################

# SAM file path
SAM_files="/scratch/groups/rosetree/RNAseq_2020/hisat2_alignment/*.sam"

# Command to convert SAM files to BAM:
for i in $SAM_files
do
tmp=${i/hisat2_alignment/BAM_files/raw_BAM}
samtools view -S -b -h $i > ${tmp/.sam/.bam}
done

# Arguments:
# -S       ignored (input format is auto-detected) - i.e. input is SAM
# -b       output in BAM format
# -h       include header in SAM output



########################################################################
#                        2. Sorting BAM file                           #
########################################################################

# BAM files need to sorted so that we can view the file in IGV 

# BAM file path - raw BAM (i.e. SAM to BAM straight conversion)
BAM_files="/scratch/groups/rosetree/RNAseq_2020/BAM_files/raw_BAM/*.bam"

# Sorting the BAM file
# Command to output a sorted file
for i in $BAM_files
do
samtools sort $i -o ${i/raw_BAM/sorted_BAM}
done

# Arguments:
# -o FILE    Write final output to FILE rather than standard output



########################################################################
#                   3. Keeping only mapped reads                       #
########################################################################

# Sorted BAM file path 
sorted_BAM_files="/scratch/groups/rosetree/RNAseq_2020/BAM_files/sorted_BAM/*.bam"

# Comand to keep on reads which have mapped/aligned to the genome
for i in $sorted_BAM_files
do
samtools view -F 4 -b $i > ${i/sorted_BAM/final_BAM}
done


# Arguments:
# -b       output BAM
# -F INT   only include reads with none of the FLAGS in INT present [0]

# Flag 4 = read unmapped (0x4)



########################################################################
#                  4. Creating indices for BAM files                   #
########################################################################
 
# We make indices so we can view the BAM files in IGV 

# Final BAM file path 
final_BAM_files="/scratch/groups/rosetree/RNAseq_2020/BAM_files/final_BAM/*.bam"

# Command to create the index file (.bai), based on the final_BAM_file 
for i in $final_BAM_files
do
samtools index $i 
done
# this creates a $final_BAM_file.bai files

# Command to creating the index file, based on the sorted_BAM_file (this one has not had reads filtered out, only sorted in chromosome position order)
for i in $sorted_BAM_files
do
samtools index $i
done
# this creates a $sorted_BAM_file.bai files



########################################################################
#                  5. FastQC of final BAM files                        #
########################################################################

# The output folder file path for the generated FastQC files
fastqc_output_folder="/scratch/groups/rosetree/RNAseq_2020/QC/FastQC_on_final_BAM_files"

# Command to perform FastQC on final BAM files
fastqc $final_BAM_files -o $fastqc_output_folder -q

# SYNOPSIS: fastqc seqfile1 seqfile2 .. seqfileN
# Arguments:
# -o --outdir     Create all output files in the specified output directory.
#                    Please note that this directory must exist as the program
#                    will not create it.  If this option is not set then the 
#                    output file for each sequence file is created in the same
#                    directory as the sequence file which was processed.
# -q --quiet      Supress all progress messages on stdout and only report errors.



########################################################################
#              Checks.... Counting reads in sorted BAM                 #
########################################################################

# Command to print total number of reads (with various restrictions detailed in the echos)
echo "\n\nPlease note that these counts are of the sorted raw BAMs before any filtering:"

for i in $sorted_BAM_files
do
echo -e "\ntotal reads of file" $i
samtools view -c $i 
echo -e "\ntotal number of properly paired reads of file" $i
samtools view -c -f 2 $i 
echo -e "\ntotal number of NOT properly paired reads of file" $i
samtools view -c -F 2 $i
echo -e "\ntotal number of mapped reads of file" $i
samtools view -c -F 4 $i
echo -e "\ntotal number of unmapped reads of file" $i
samtools view -c -f 4 $i 
echo -e "\ntotal number of reads that pass quality score 20 of file" $i
samtools view -c -q 20 $i 
echo -e "\n\n"
done

# Arguments:
# -f INT   only include reads with all  of the FLAGs in INT present [0]
# -F INT   only include reads with none of the FLAGS in INT present [0]
# -q INT   only include reads with mapping quality >= INT [0]
# -c       print only the count of matching records

# If you just type "samtools view file.bam" w/out the -c it prints the entire file to your screen!! 



########################################################################
#              Checks.... Verifying sorting was successful             #
########################################################################
 
# Previewing the sorted BAMs

# Command to print BAM file head & tail
for i in $final_BAM_files
do
echo -e "\n\nBAM file head of" $i 
samtools view $i | head -n 3
echo -e "\n\nBAM file tail of" $i 
samtools view $i | tail -n 3
echo -e "\n\n\n"
done



########################################################################
#                             The end                                  #
########################################################################

# print end time to slurm output file
echo -e "\nend: " `date`
