# README

## Description
This repository contains supporting materials for an upcoming article. It includes R and Python code for reproducing the results presented in the article.

## Reproducing the Results
To reproduce the results obtained from the R code, follow these steps:

### Install R
If you haven't already, you'll need to install R. You can download it from the [official R website](https://www.r-project.org/).

### 1. Linkage disequilibrium Pruning
#### 1.1 Install Required Packages
Make sure you have the required R packages installed. You can install them by running the following command in your R console:

```R
install.packages(c("pedtools", "devtools", "usethis"))
devtools::install_github("magnusdv/pedbuildr")
devtools::install_github("magnusdv/forrel")
```

#### 1.2 Input files
This script utilizes two input files:
1. A file containing genetic information of individuals, following the schema below:

| id | fid | mid | sex | rs1 | rs2 | rs3 | rs4 | rs5 |
|----|-----|-----|-----|-------|-------|-------|-------|-------|
| 1  | 0   | 0   | 1   | 1/2   | 1/2   | 1/1   | 1/2   | 1/2   |
| 2  | 0   | 0   | 2   | 1/1   | 1/2   | 1/2   | 2/2   | 1/2   |
| 3  | 1   | 2   | 1   | 1/1   | 1/2   | 1/2   | 2/2   | 1/1   |


- `id`: Individual ID
- `fid`: Father's ID
- `mid`: Mother's ID
- `sex`: Sex (1 for male, 2 for female)
- `<rs1>`, `<rs2>`, `<rs3>`, etc.: Genotypes for each marker (e.g., 1/2, 1/1, 2/2)

2. A file containing information about markers (RS IDs and allele frequencies), following the schema below:

alleles|frequencies
----|----
rs1 | NA
1   | 0.5
2   | 0.5
rs2 | NA
1   | 0.75
2   | 0.25
rs3 | NA
1   | 0.33
2   | 0.33
3   | 0.33
  

#### 2. IBD and LR estimation
#### 2.1 Install Required Packages
Make sure you have the required R packages installed. You can install them by running the following command in your R console:
```R
install.packages(c("LDlinkR", "remotes"))
remotes::install_github("CBIIT/LDlinkR")
```
