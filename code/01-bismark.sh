#!/bin/bash
## Job Name
#SBATCH --job-name=Cvbs
## Allocation Definition
#SBATCH --account=coenv
#SBATCH --partition=coenv
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=10-00:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr320@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/sr320/120321-cvBS



# Directories and programs
bismark_dir="/gscratch/srlab/programs/Bismark-0.22.3"
bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/"
samtools="/gscratch/srlab/programs/samtools-1.9/samtools"
reads_dir="/gscratch/srlab/sr320/data/cvirg/"
fastqc="/gscratch/srlab/programs/fastqc_v0.11.9/fastqc"
genome_folder="/gscratch/srlab/sr320/data/Cvirg-genome/"

source /gscratch/srlab/programs/scripts/paths.sh

# 
# 
# ${bismark_dir}/bismark_genome_preparation \
# --verbose \
# --parallel 28 \
# --path_to_aligner ${bowtie2_dir} \
# ${genome_folder}


# 
# FastQC 
# 
# Populate array with FastQ files
#fastq_array=(/gscratch/srlab/sr320/data/cvirg/*.fq.gz)
# 
# Pass array contents to new variable
#fastqc_list=$(echo "${fastq_array[*]}")
# 
# Run FastQC
#${fastqc} \
#--threads 28 \
#${fastqc_list} \
#-o /gscratch/scrubbed/sr320/071521-cvBS/




# zr3534_10_R1.fq.gz  zr3644_11_R2.fq.gz  zr3644_22_R1.fq.gz
# 0502_R1.fq.gz  zr3534_10_R2.fq.gz  zr3644_12_R1.fq.gz  zr3644_22_R2.fq.gz

# find ${reads_dir}*_1.fq.gz \
# | xargs basename -s _R1_val_1.fq.gz | xargs -I{} ${bismark_dir}/bismark \
# --path_to_bowtie ${bowtie2_dir} \
# -genome ${genome_folder} \
# -p 4 \
# -score_min L,0,-0.6 \
# --non_directional \
# -1 ${reads_dir}{}_R1_val_1.fq.gz \
# -2 ${reads_dir}{}_R2_val_2.fq.gz
# 
# find *.bam | \
# xargs basename -s .bam | \
# xargs -I{} ${bismark_dir}/deduplicate_bismark \
# --bam \
# --paired \
# {}.bam


# 
# ${bismark_dir}/bismark_methylation_extractor \
# --bedGraph \
# --counts \
# --comprehensive \
# --merge_non_CpG \
# --multicore 28 \
# --buffer_size 75% \
# *deduplicated.bam


# 
# # Bismark processing report
# 
# ${bismark_dir}/bismark2report
# 
# #Bismark summary report
# 
# ${bismark_dir}/bismark2summary
# 
#  #run multiqc
# /gscratch/srlab/programs/anaconda3/bin/multiqc .
# 

# Sort files for methylkit and IGV

find *deduplicated.bam | \
xargs basename -s .bam | \
xargs -I{} ${samtools} \
sort --threads 28 {}.bam \
-o {}.sorted.bam

# Index sorted files for IGV
# The "-@ 16" below specifies number of CPU threads to use.

find *.sorted.bam | \
xargs basename -s .sorted.bam | \
xargs -I{} ${samtools} \
index -@ 28 {}.sorted.bam





find *deduplicated.bismark.cov.gz \
| xargs basename -s _bismark_bt2_pe.deduplicated.bismark.cov.gz \
| xargs -I{} ${bismark_dir}/coverage2cytosine \
--genome_folder ${genome_folder} \
-o {} \
--merge_CpG \
--zero_based \
{}_bismark_bt2_pe.deduplicated.bismark.cov.gz


#creating bedgraphs post merge

for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4}}' \
  > "${STEM}"_10x.bedgraph
done



for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4}}' \
  > "${STEM}"_5x.bedgraph
done


#creating tab files with raw count for glms

for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_10x.tab
done


for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_5x.tab
done


