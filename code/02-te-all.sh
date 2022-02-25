#!/bin/bash
## Job Name
#SBATCH --job-name=Cov
## Allocation Definition
#SBATCH --account=coenv
#SBATCH --partition=coenv
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=5-00:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr320@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/sr320/022122-te-cov



# Directories and programs
bismark_dir="/gscratch/srlab/programs/Bismark-0.22.3"
bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/"
samtools="/gscratch/srlab/programs/samtools-1.9/samtools"
reads_dir="/gscratch/srlab/sr320/data/cvirg/"
fastqc="/gscratch/srlab/programs/fastqc_v0.11.9/fastqc"
genome_folder="/gscratch/srlab/sr320/data/Cvirg-genome/"

source /gscratch/srlab/programs/scripts/paths.sh

cd /gscratch/scrubbed/sr320/021022-CvBSnp/
FILES=$(ls *sorted.bam)
cd /gscratch/scrubbed/sr320/022122-te-cov/

for file in ${FILES}
do
    NAME=$(echo ${file} | awk -F "." '{print $1}')
    echo ${NAME}

    /gscratch/srlab/programs/bedtools-2.27.1/bin/bedtools coverage \
    -mean \
    -sorted \
    -a /gscratch/srlab/sr320/data/Cvirg-genome/C_virginica-3.0_TE-all.gff \
    -b /gscratch/scrubbed/sr320/021022-CvBSnp/${NAME}.deduplicated.sorted.bam \
    > ${NAME}_te-allcov_m.out
done
