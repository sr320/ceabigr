# 18-differential-gene-expression-ballgown

`ceabigr/output/18-differential-gene-expression-ballgown`

---

All files in this directory generated by [`18-differential-gene-expression-ballgown.Rmd`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown).

- [`DEG_sex_female_filtered_p-0.05_q-0.05.bed`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DEG_sex_female_filtered_p-0.05_q-0.05.bed): BED file of differentially expressed transcripts (DEGs) up-regulated in females only, filtered for DEGs having p-values and q-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.

- [`DEG_sex_filtered_p-0.05_q-0.05.bed`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DEG_sex_filtered_p-0.05_q-0.05.bed): BED file of differentially expressed genes (DEGs) identified between sexes, filtered for DEGs having p-values and q-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.

- [`DEG_sex_filtered_p-0.05_q-0.05.csv`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DEG_sex_filtered_p-0.05_q-0.05.csv): Differentially expressed genes (DEGs) identified between sexes, filtered for DEGs having p-values and q-values <= 0.05. Contains the following columns:

  - `t_id`: Transcript ID (integer) assigned by `ballgown`.

  - `chr`: Chromosome ID from NCBI [_Crassostrea virginica_ (Eastern oyster)](https://en.wikipedia.org/wiki/Eastern_oyster) genome `GCF_002022765.2_C_virginica-3.0_genomic.fna`.

  - `strand`: Originating DNA strand (+ or -).

  - `start`: One-based starting location of transcript on chromosome.

  - `end`: One-based ending location of transcript on chromosome.

  - `t_name`: Annotated transcript name taken from genome file.

  - `num_exons`: Number of exons comprising transcript.

  - `length`: Length of transcript (bp).

  - `gene_id`: Corresponding gene ID, if applicable.

  - `gene_name`: Corresponding gene name, if applicable.

  - `cov.*`: Sequencing read coverage for transcript for each sample.

  - `FPKM.*`: Relative expression levels for each transcript in each sample.

- [`DEG_sex_male_filtered_p-0.05_q-0.05.bed`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DEG_sex_male_filtered_p-0.05_q-0.05.bed): BED file of differentially expressed genes (DEGs) up-regulated in males only, filtered for DEGs having p-values and q-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.

- [`DET_females_controls_p-0.05_q-0.05.bed`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DET_females_controls_p-0.05_q-0.05.bed): BED file of differentially expressed transcripts (DETs) up-regulated in female control samples only, filtered for DETs having p-values and q-values <= 0.05. Includes optional columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.

- [`DET_females_controls_vs_exposed_p-0.05_q-0.05.bed`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DET_females_controls_vs_exposed_p-0.05_q-0.05.bed): BED file of differentially expressed transcripts (DETs) identified between females controls and exposed samples, filtered for DETs having p-values and q-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.

- [`DET_females_controls_vs_exposed_p-0.05_q-0.05.csv`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DET_females_controls_vs_exposed_p-0.05_q-0.05.csv): Differentially expressed transcripts (DETs) identified between females controls and exposed samples, filtered for DETs having p-values and q-values <= 0.05. Contains the following columns:

  - `t_id`: Transcript ID (integer) assigned by `ballgown`.

  - `chr`: Chromosome ID from NCBI [_Crassostrea virginica_ (Eastern oyster)](https://en.wikipedia.org/wiki/Eastern_oyster) genome `GCF_002022765.2_C_virginica-3.0_genomic.fna`.

  - `strand`: Originating DNA strand (+ or -).

  - `start`: One-based starting location of transcript on chromosome.

  - `end`: One-based ending location of transcript on chromosome.

  - `t_name`: Annotated transcript name taken from genome file.

  - `num_exons`: Number of exons comprising transcript.

  - `length`: Length of transcript (bp).

  - `gene_id`: Corresponding gene ID, if applicable.

  - `gene_name`: Corresponding gene name, if applicable.

  - `cov.*`: Sequencing read coverage for transcript for each sample.

  - `FPKM.*`: Relative expression levels for each transcript in each sample.

- [`DET_females_exposed_p-0.05_q-0.05.bed`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DET_females_exposed_p-0.05_q-0.05.bed): BED file of differentially expressed transcripts (DETs) up-regulated in females exposed samples, filtered for DETs having p-values and q-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.

- [`DET_sex_female_filtered_p-0.05_q-0.05.bed`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DET_sex_female_filtered_p-0.05_q-0.05.bed): BED file of differentially expressed transcripts (DETs) up-regulated in females only, filtered for DETs having p-values and q-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.

- [`DET_sex_filtered_p-0.05_q-0.05.bed`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DET_sex_filtered_p-0.05_q-0.05.bed): BED file of differentially expressed transcripts (DETs) identified between sexes, filtered for DETs having p-values and q-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.

- [`DET_sex_filtered_p-0.05_q-0.05.csv`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DET_sex_filtered_p-0.05_q-0.05.csv): Differentially expressed transcripts (DETs) identified between sexes, filtered for DETs having p-values and q-values <= 0.05. Contains the following columns:

  - `t_id`: Transcript ID (integer) assigned by `ballgown`.

  - `chr`: Chromosome ID from NCBI [_Crassostrea virginica_ (Eastern oyster)](https://en.wikipedia.org/wiki/Eastern_oyster) genome `GCF_002022765.2_C_virginica-3.0_genomic.fna`.

  - `strand`: Originating DNA strand (+ or -).

  - `start`: One-based starting location of transcript on chromosome.

  - `end`: One-based ending location of transcript on chromosome.

  - `t_name`: Annotated transcript name taken from genome file.

  - `num_exons`: Number of exons comprising transcript.

  - `length`: Length of transcript (bp).

  - `gene_id`: Corresponding gene ID, if applicable.

  - `gene_name`: Corresponding gene name, if applicable.

  - `cov.*`: Sequencing read coverage for transcript for each sample.

  - `FPKM.*`: Relative expression levels for each transcript in each sample.

- [`DET_sex_male_filtered_p-0.05_q-0.05.bed`](https://github.com/sr320/ceabigr/tree/main/output/18-differential-gene-expression-ballgown/DET_sex_male_filtered_p-0.05_q-0.05.bed): BED file of differentially expressed transcripts (DETs) up-regulated in males only, filtered for DETs having p-values and q-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.