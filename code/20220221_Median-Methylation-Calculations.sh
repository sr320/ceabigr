#!/bin/bash
## Job Name
#SBATCH --job-name=BmrkAln_Salmo
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=5-23:30:00
## Memory per node
#SBATCH --mem=500G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/strigg/analyses/20200613


%%bash

#record each command that is called in slurm file

set -ex

#add paths to programs to bash profile

source /gscratch/srlab/strigg/bin/paths.sh

#print the path of the programs

which samtools
which bowtie2
which bismark

#define genome folder

genome_folder="/gscratch/srlab/strigg/data/Ssalar/GENOMES/chr1-29MT"

#run bismark in PE mode on trimmed reads

find /gscratch/scrubbed/strigg/analyses/20200612/*_R1_001_val_1.fq.gz |\
xargs basename -s _R1_001_val_1.fq.gz| \
xargs -I{} bismark \
--score_min L,0,-0.3 \
-p 4 \
--non_directional \
--genome ${genome_folder} \
-1 /gscratch/scrubbed/strigg/analyses/20200612/{}_R1_001_val_1.fq.gz \
-2 /gscratch/scrubbed/strigg/analyses/20200612/{}_R2_001_val_2.fq.gz \

# deduplicate

find *.bam | \
xargs basename -s .bam | \
xargs -I{} deduplicate_bismark \
--bam \
--paired \
{}.bam

#run methylation extractor
bismark_methylation_extractor \
--paired-end \
--bedGraph \
--comprehensive \
--counts \
--scaffolds \
--multicore 14 \
--buffer_size 75% \
*.deduplicated.bam


#######################################################
### Generate bed graphs for IGV and bedtools analysis##
#######################################################

# Sort files for methylkit and IGV

find *.deduplicated.bam | \
xargs basename -s .bam | \
xargs -I{} samtools \
sort --threads 28 {}.bam \
-o {}.sorted.bam

# Index sorted files for IGV
# The "-@ 16" below specifies number of CPU threads to use.

find *.sorted.bam | \
xargs basename -s .sorted.bam | \
xargs -I{} samtools \
index -@ 28 {}.sorted.bam

# Run multiqc

/gscratch/srlab/strigg/bin/anaconda3/bin/multiqc \
/gscratch/scrubbed/strigg/analyses/20200612/ \
/gscratch/scrubbed/strigg/analyses/20200613/

# create merged
find *.deduplicated.bismark.cov.gz |\
xargs basename -s .deduplicated.bismark.cov.gz |\
xargs -I{} coverage2cytosine \
--genome_folder ${genome_folder} \
-o {} \
--merge_CpG \
--zero_based \
{}.deduplicated.bismark.cov.gz


#creating tab files with % me, raw mCpG and total CpG counts

for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4, $5, $5+$6}}' \
  > "${STEM}"_5x.bed
done
