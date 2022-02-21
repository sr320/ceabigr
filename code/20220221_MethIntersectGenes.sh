#!/bin/bash
## Job Name
#SBATCH --job-name=GeneMethClass
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
#SBATCH --mail-user=shellywanamaker@gmail.com
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/strigg/analyses/20220221


%%bash

#record each command that is called in slurm file

set -ex

#add paths to programs to bash profile

source /gscratch/srlab/strigg/bin/paths.sh

#define genome folder

genome_folder="/gscratch/srlab/strigg/data/Ssalar/GENOMES/chr1-29MT"

#get data
wget -r \
--no-check-certificate \
--quiet \
--no-directories --no-parent \
-P /gscratch/scrubbed/strigg/data/20220221 \  #where files go - indicating current directory
-A https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/120321-cvBS/*10x.bedgraph  \ #which files you want


#get gene.bed file
wget -r \
--no-check-certificate \
--quiet \
--no-directories --no-parent \
-P /gscratch/scrubbed/strigg/data/20220221 \  #where files go - indicating current directory
-A https://eagle.fish.washington.edu/Cvirg_tracks/C_virginica-3.0_Gnomon_genes.bed  \ #which files you want


#For each sample-gene overlap file (ends in *bedgraph.bed-mcGenes)
#Use intersectBed to find where loci and genes intersect, allowing loci to be mapped to genes
#wb: Print all lines in the second file
#a: sample-gene overlap file
#b: gene track
#Save output in a new file that has the same base name and ends in -geneID
for f in /gscratch/scrubbed/strigg/data/20220221/*10x.bedgraph
do
  intersectBed \
  -wb \
  -a ${f} \
  -b /gscratch/scrubbed/strigg/data/20220221/C_virginica-3.0_Gnomon_genes.bed \
  > ${f}-geneID
done
