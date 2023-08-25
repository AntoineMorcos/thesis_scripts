#!/bin/bash -l

#SBATCH --output=slurm-%u-%j-featureCounts.out  ## output file name
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
module add utilities/use.dev
module load apps/subread/2.0.0 



########################################################################
#                              File paths                              #
########################################################################

# File path of final BAM files
final_BAM_files="/scratch/groups/rosetree/RNAseq_2020/BAM_files/final_BAM/*.bam"

# File path of output file including read counts
output_TXT_file="/scratch/groups/rosetree/RNAseq_2020/gene_counts/counts_paired_hg38.txt"

# File path of GTF file
GTF_file="/scratch/groups/rosetree/RNAseq_2020/GTF/hg38.ncbiRefSeq.gtf" 



########################################################################
#                   Counting the reads over genes                      #
########################################################################

# Command to count reads over the genes (transcript quantification) using featureCounts
featureCounts -p -B -t exon -g gene_id -F GTF -a $GTF_file -o $output_TXT_file $final_BAM_files


# featureCounts (Version 2.0.0)

# Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ..
#  input_file1 [input_file2] ...   A list of SAM or BAM format files. They can be
#                      either name or location sorted. If no files provided,
#                      <stdin> input is expected. Location-sorted paired-end reads
#                      are automatically sorted by read names.

# Arguments:
#  -p                  If specified, fragments (or templates) will be counted
#                      instead of reads. This option is only applicable for
#                      paired-end reads; single-end reads are always counted as
#                      reads.
#
#  -B                  Only count read pairs that have both ends aligned.
#
#  -t <string>         Specify feature type in GTF annotation. 'exon' by 
#                      default. Features used for read counting will be 
#                      extracted from annotation using the provided value.
#
#  -g <string>         Specify attribute type in GTF annotation. 'gene_id' by 
#                      default. Meta-features used for read counting will be 
#                      extracted from annotation using the provided value.
#
#  -a <string>         Name of an annotation file. GTF/GFF format by default. See
#                      -F option for more format information. Inbuilt annotations
#                      (SAF format) is available in 'annotation' directory of the
#                      package. Gzipped file is also accepted.
#
#  -F <string>         Specify format of the provided annotation file. Acceptable
#                      formats include 'GTF' (or compatible GFF format) and
#                      'SAF'. 'GTF' by default.  For SAF format, please refer to
#                      Users Guide.
#                      
#  -o <string>         Name of output file including read counts. A separate file
#                      including summary statistics of counting results is also
#                      included in the output ('<string>.summary'). Both files
#                      are in tab delimited format.



########################################################################
#                             The end                                  #
########################################################################

# print end time to slurm output file
echo -e "\nend: " `date`
