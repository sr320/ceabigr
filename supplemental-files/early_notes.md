# NOTES

## HYPOTHESES

- **Epigenetics (DNA methylation) modulates genome to phenome linkages**
- **Stressor perturbs the system**

### Set the stage/explain the -omic landscape (briefly):
- Genetic relatedness
- DEG
- DMG

### “Assign” variation to genetic or epigenetic
- Epigenetic WRT genetic influence hypotheses

## Genetic Influence (TE, SNPs, not limited to same gene)
- Lower transposable element activity with higher methylation activity
- Individual genetic differences are associated with methylation differences
- SNPs associated with methylation change (QTL)
- CpG-SNPs (create or remove a CpG, change the possibility of methylation occurring)
- Next step: Get list of methylated sites that have a genetic influence
- There are differentially methylated loci influenced by the environment that are not associated with genetic differences
- Next step: take mQTL from (2) and removing them from DML for sex and treatment
- Methylation leads to genetic variation that is beneficial
- More methylation → more CpGs → less genetic variation outside of CpGs
- Can determine if we see patterns of genetic variation/methylation

### Role of DNA Methylation on Gene Activity (expression, splice variants, lncRNA, copy number)
- There is a positive relationship between gene expression and methylation
- Methylation modulates phenotype via gene expression
- Gene expression ~ methylation at each exon (data available, in Javie and Ariana’s Raven repositories)
- Next step: Subset/investigate in subsets of categories/families of genes and run regressions/correlations for each
- Next step: SPLS analysis with mixOmics package
- Next step: Run these analyses in DEGs/DMGs
- There is negative relationship between DNA methylation and CoV of gene expression - DONE with transcriptional noise analysis
- We have gene expression data and methylation for each exon (data available, in Javie and Ariana’s Raven repositories)
- Next step: Calculate CoV of each exon/gene unit of interest and run analyses
- Next step: Run these analyses in DEGs/DMGs
- There is a negative relationship between methylation status and alternative splicing
- \# transcripts as measure for alternative splicing 
- Next step: run correlations between transcript and isoform count information
- Next step: Run these analyses in DEGs/DMGs
- Exon specific methylation influences quantitative expression or specific isoforms
- Long non-coding RNAs and methylation

### Synthesis: Connections between genetics, DNA Methylation, and Gene Activity (eQTLs, multi-omic correlations/networks, QTL between methylation and expression, TE-> methylation -> expression )

#### Approaches:
- Intersecting gene sets from various analyses
- QTL analyses:  Identify methylated loci, SNPs, and expressed genes that are all commonly associated locally and distally. 
- Mixomics package 
- Other multi-omic package that considers genotype 

- Separation of groups visualized on block PLSDA - how do individual -omics separate groups compared to the multi-omic analysis? 
- What are the primary drivers of the separation? E.g., are differences by sex driven by methylation or are differences by treatment driven by gene expression?
- Are there correlations between methylation and gene expression (e.g., view on a circos plot) 
- There will be differences in phenotype in response to environment between male and female gametes

---

## Overview
### Genetic relatedness
- Todo: generate high quality SNP file, generate genetic relatedness matrix
### DEG
### DMG - 
- Hollie Putnam was working on the DMG analysis for egg samples, sperm samples, and the whole dataset. These are calculated for CpGs with 10x coverage and are found in all samples in the comparison being made.  These filtered CpGs were also used to calculate the mean methylation per gene.
- ToDo finalize these 3 outputs and figures

## Genetic Influence (TE, SNPs, not limited to same gene)
- Lower transposable element activity with higher methylation activity
- Generated single base depth files from WGBS .bam files with samtools, 10x methylation .cov files and expression per base depth with bedtools. Methylation level per TE type was calculated and a graph was generated per type and treatment. Expression and copy number analyses are still missing (running). 
 
- Individual genetic differences are associated with methylation differences
- KES: started methQTL analysis with matrixEQTL in R
- Todo: identify methQTLs associated with DMLs (by treatment, indv, sex)
- Todo: identify CpG-SNPs and associate with DMLs, methQTLs
- Todo: correlate
