`ceabigr/output/17.1-GO-and-GOslim-CDS-Annotation`

Output files generated by [`17.1-GO-and-GOslim-CDS-Annotation.Rmd`](https://github.com/sr320/ceabigr/blob/main/code/17.1-GO-and-GOslim-CDS-Annotation.Rmd).

- [`Cvir-CDS-geneID-GOID.tab`](https://github.com/sr320/ceabigr/tree/main/output/17.1-GO-and-GOslim-CDS-Annotation/Cvir-CDS-geneID-GOID.tab): Two column, tab-delimted file of genes and corresponding GO IDs.
GO IDs are separated by "; ". Note the space after the semi-colon. 

- [`Cvir-CDS-GOID-geneID.tab`](https://github.com/sr320/ceabigr/tree/main/output/17.1-GO-and-GOslim-CDS-Annotation/Cvir-CDS-GOID-geneID.tab): Two column, tab-delimted file of GO IDs and corresponding gene.
GO IDs are separated by ";".

- [`Cvir-CDS-GOslim.BP_per_gene.tab`](https://github.com/sr320/ceabigr/tree/main/output/17.1-GO-and-GOslim-CDS-Annotation/Cvir-CDS-GOslim.BP_per_gene.tab): Three column, tab-delimited file of genes and
corresponding Biological Process GOslim(s) and GOslim term(s). GOslims and GOslim terms are semi-colon delimited in instances where there are more than one for a given gene.

- [`Cvir-CDS-GOslim.BP_term_GOIDs_genes.tab`](https://github.com/sr320/ceabigr/tree/main/output/17.1-GO-and-GOslim-CDS-Annotation/Cvir-CDS-GOslim.BP_term_GOIDs_genes.tab): Four column, tab-delimited file with one
Biological Process GOslim ID per row. GO IDs are semi-colon delimited. Genes are comma-delimited.
Columns are: GOslim.BP, Term, GO.IDs, Genes

- [`Cvir-CDS-uniprot-full.tsv`](https://github.com/sr320/ceabigr/tree/main/output/17.1-GO-and-GOslim-CDS-Annotation/Cvir-CDS-uniprot-full.tsv): Tab-delimited output file from UniProt API retrieval.
Columns:
  - `Entry`
  - `Reviewed`
  - `Entry Name`
  - `Protein names`
  - `Gene Names`
  - `Organism`
  - `Length`
  - `Gene Ontology (biological process)`
  - `Gene Ontology (cellular component)`
  - `Gene Ontology (GO)`
  - `Gene Ontology (molecular function)`
  - `Gene Ontology IDs`
