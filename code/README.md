# code

Code used in for each analysis can be found below.

## Genome Information and Feature Tracks

- [`Generating-Genome-Feature-Tracks.ipynb`](https://github.com/sr320/ceabigr/blob/main/code/Generating-Genome-Feature-Tracks.ipynb): Generate genome feature tracks based on NCBI annotations. Output can be fond in the [`genome-features`](https://github.com/sr320/ceabigr/tree/main/genome-features) directory

## Genetic Variation

- [`20220921-cvir-ceabigr-nextflow-epidiverse-snp.sh`](https://github.com/sr320/ceabigr/blob/main/code/20220921-cvir-ceabigr-nextflow-epidiverse-snp.sh): Bash script to run the Nextflow EpiDiverse/snp pipeline.

## Gene Activity

### Gene and Transcript Expression

- [`20220131_cvir_hisat2-GCF_002022765.2_adult-oa-gonad.sh`](https://github.com/sr320/ceabigr/blob/main/code/20220131_cvir_hisat2-GCF_002022765.2_adult-oa-gonad.sh): Bash script to align RNA-seq data to _C.virginica_ genome using HISAT2.
- [`18-differential-gene-expression-ballgown.Rmd`](https://github.com/sr320/ceabigr/blob/main/code/18-differential-gene-expression-ballgown.Rmd): Differential gene and transcript expression analysis using `ballgown` R package.
- [`00-75-max-transcript-enrichment.qmd`](https://github.com/sr320/ceabigr/blob/main/code/00-75-max-transcript-enrichment.qmd): Identify genes with changes in the maximum number of transcripts expressed

### Predominant Transcript Identification

- [`42-predominant-isoform.Rmd)`](https://github.com/sr320/ceabigr/blob/main/code/00-42-predominant-isoform.Rmd): Identify genes with a shift in the predominant transcript

### Alternative Splicing

## DNA Methylation Analysis

- [`bismark`](https://github.com/sr320/ceabigr/blob/main/code/01-bismark.sh): bismark code
- [`General-Methylation-Landscape.ipynb`](https://github.com/sr320/ceabigr/blob/main/code/General-Methylation-Landscape.ipynb): Use `bedtools` to characterize methylation landscape using CpGs with 10x coverage
- [`General-Methylation-Landscape.Rmd`](https://github.com/sr320/ceabigr/blob/main/code/General-Methylation-Landscape.Rmd): Visualization of general methylation landscape and chi-squared tests

### CpG locus-level Methylation Quantification

- [`Genomic-Location-of-DML.ipynb`](https://github.com/sr320/ceabigr/blob/main/code/00-Genomic-Location-of-DML.ipynb): Use `bedtools` to identify genomic locations of DML
- [`Genomic-Location-of-DML.Rmd`](https://github.com/sr320/ceabigr/blob/main/code/Genomic-Location-of-DML.Rmd): Visualization of DML genomic locations and chi-squared tests

### Gene-level Methylation Quantification

## Methylation Influence on Gene Activity

- [`00-75-max-transcript-enrichment.qmd`](https://github.com/sr320/ceabigr/blob/main/code/00-75-max-transcript-enrichment.qmd): Understand influence of methylation on genes with changes in the maximum number of transcripts expressed and conduct enrichment analyses

## Influence of Methylation on Predominant Transcript Shifts

- [`42-predominant-isoform.Rmd)`](https://github.com/sr320/ceabigr/blob/main/code/00-42-predominant-isoform.Rmd): Understand the influence of methylation on predominant transcript shifts and conduct enrichment analyses

## Transcriptional Noise Analysis

- [`42-predominant-isoform.Rmd)`](https://github.com/sr320/ceabigr/blob/main/code/00-42-predominant-isoform.Rmd): Transcriptional noise analysis

## Influence of Methylation on Alternative Splicing

- [`00-78-asca-methdiff.qmd`](https://github.com/sr320/ceabigr/blob/main/code/00-78-asca-methdiff.qmd): Understand if methylation significantly impacted alternative splicing and conduct enrichment analyses

## Other scripts

- [`01.01-supplemental-tables.Rmd`](01.01-supplemental-tables.Rmd): R Markdown document for generating supplemental tables for manuscript. Output files are in [`supplemental-files/`](../supplemental-files/)
- [`98-enrichment-result-compilation.Rmd`](https://github.com/sr320/ceabigr/blob/main/code/98-enrichment-result-compilation.Rmd): Compilation of enrichment results across all tests
