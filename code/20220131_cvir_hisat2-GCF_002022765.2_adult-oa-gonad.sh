#!/bin/bash
## Job Name
#SBATCH --job-name=20220131_cvir_hisat2-GCF_002022765.2_adult-oa-gonad
## Allocation Definition
#SBATCH --account=coenv
#SBATCH --partition=coenv
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=5-00:00:00
## Memory per node
#SBATCH --mem=200G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samwhite@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/samwhite/outputs/20220131_cvir_hisat2-GCF_002022765.2_adult-oa-gonad

## Hisat2 alignment of C.virginica adult OA gonad RNAseq to NCBI C.virginica genome assembly
## using HiSat2 index generated on 20210720.

## Expects FastQ input filenames to match *fastp-trim.20210714.fq.gz


###################################################################################
# These variables need to be set by user

## Assign Variables

# Set number of CPUs to use
threads=40

# Index name for Hisat2 use
# Needs to match index naem used in previous Hisat2 indexing step
genome_index_name="cvir_GCF_002022765.2"

# Set output filename
sample_name="20210131-cvir-hisat2"

# Paths to programs
hisat2_dir="/gscratch/srlab/programs/hisat2-2.1.0"
hisat2="${hisat2_dir}/hisat2"
samtools="/gscratch/srlab/programs/samtools-1.10/samtools"


# Input/output files
genome_index_dir="/gscratch/srlab/sam/data/C_virginica/genomes"
fastq_dir="/gscratch/srlab/sam/data/C_virginica/RNAseq/"


# Programs associative array
declare -A programs_array
programs_array=(
[hisat2]="${hisat2}" \
[samtools_index]="${samtools} index" \
[samtools_sort]="${samtools} sort" \
[samtools_view]="${samtools} view"
)


###################################################################################################

# Exit script if any command fails
set -e

# Load Python Mox module for Python module availability

module load intel-python3_2017

## Inititalize arrays
fastq_array_R1=()
fastq_array_R2=()

# Copy Hisat2 genome index files
rsync -av "${genome_index_dir}"/${genome_index_name}*.ht2 .

# Create array of fastq R1 files
# and generate MD5 checksums file.
for fastq in "${fastq_dir}"*fastp-trim.20210714.fq.gz
do
  fastq_array_R1+=("${fastq}")
  echo "Generating checksum for ${fastq}..."
  md5sum "${fastq}" >> input_fastqs_checksums.md5
  echo "Checksum for ${fastq} completed."
  echo ""
done

# Create array of fastq R2 files
# and append to MD5 checksums file.
for fastq in "${fastq_dir}"*fastp-trim.20210714.fq.gz
do
  fastq_array_R2+=("${fastq}")
  echo "Generating checksum for ${fastq}..."
  md5sum "${fastq}" >> input_fastqs_checksums.md5
  echo "Checksum for ${fastq} completed."
  echo ""
done

# Create comma-separated lists of FastQs for Hisat2
printf -v joined_R1 '%s,' "${fastq_array_R1[@]}"
fastq_list_R1=$(echo "${joined_R1%,}")

printf -v joined_R2 '%s,' "${fastq_array_R2[@]}"
fastq_list_R2=$(echo "${joined_R2%,}")


# Hisat2 alignments
"${programs_array[hisat2]}" \
-x "${genome_index_name}" \
-1 "${fastq_list_R1}" \
-2 "${fastq_list_R2}" \
-S "${sample_name}".sam \
2> "${sample_name}"_hisat2.err

# Sort SAM files, convert to BAM, and index
${programs_array[samtools_view]} \
-@ "${threads}" \
-Su "${sample_name}".sam \
| ${programs_array[samtools_sort]} - \
-@ "${threads}" \
-o "${sample_name}".sorted.bam
${programs_array[samtools_index]} "${sample_name}".sorted.bam


# Delete unneccessary index files
rm "${genome_index_name}"*.ht2

# Delete unneded SAM files
rm ./*.sam

# Generate checksums
for file in *
do
  md5sum "${file}" >> checksums.md5
done

#######################################################################################################

# Capture program options
if [[ "${#programs_array[@]}" -gt 0 ]]; then
  echo "Logging program options..."
  for program in "${!programs_array[@]}"
  do
    {
    echo "Program options for ${program}: "
    echo ""
    # Handle samtools help menus
    if [[ "${program}" == "samtools_index" ]] \
    || [[ "${program}" == "samtools_sort" ]] \
    || [[ "${program}" == "samtools_view" ]]
    then
      ${programs_array[$program]}

    # Handle DIAMOND BLAST menu
    elif [[ "${program}" == "diamond" ]]; then
      ${programs_array[$program]} help

    # Handle NCBI BLASTx menu
    elif [[ "${program}" == "blastx" ]]; then
      ${programs_array[$program]} -help
    fi
    ${programs_array[$program]} -h
    echo ""
    echo ""
    echo "----------------------------------------------"
    echo ""
    echo ""
  } &>> program_options.log || true

    # If MultiQC is in programs_array, copy the config file to this directory.
    if [[ "${program}" == "multiqc" ]]; then
      cp --preserve ~/.multiqc_config.yaml multiqc_config.yaml
    fi
  done
fi


# Document programs in PATH (primarily for program version ID)
{
date
echo ""
echo "System PATH for $SLURM_JOB_ID"
echo ""
printf "%0.s-" {1..10}
echo "${PATH}" | tr : \\n
} >> system_path.log