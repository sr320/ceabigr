#!/bin/bash
## Job Name
#SBATCH --job-name=haplo
## Allocation Definition
#SBATCH --account=coenv
#SBATCH --partition=coenv
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=20-00:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr320@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/sr320/031722-dbl-thread



# Directories and programs
bismark_dir="/gscratch/srlab/programs/Bismark-0.22.3"
bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/"
samtools="/gscratch/srlab/programs/samtools-1.9/samtools"
reads_dir="/gscratch/srlab/sr320/data/cvirg/"
fastqc="/gscratch/srlab/programs/fastqc_v0.11.9/fastqc"
genome_folder="/gscratch/srlab/sr320/data/Cvirg-genome/"

source /gscratch/srlab/programs/scripts/paths.sh

for file in ../031722-haplo/*dedup-split-RG.bam
do
sample="$(basename -a $file | cut -d "." -f 1)"

/gscratch/srlab/programs/gatk-4.1.9.0/gatk HaplotypeCaller \
-R ../031722-haplo/Cvirginica_v300.fa \
-I ../031722-haplo/$sample.dedup-split-RG.bam \
-O $sample.variants.g.vcf \
--native-pair-hmm-threads 28 \
-ERC GVCF
done

