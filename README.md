# README

## Description

This repository contains supporting materials for an upcoming article. It includes R and Python code for reproducing the results presented in the article.

## R 📊 scripts

### Install R

If you haven't already, you'll need to install R. You can download it from the [official R website](https://www.r-project.org/).

---

### 1. Linkage disequilibrium Pruning 🔗✂️

This script is designed to identify independent genetic variants by assessing linkage disequilibrium.

#### 1.1 Install Required Packages

Make sure you have the required R packages installed. You can install them by running the following command in your R console:

```R
install.packages(c("LDlinkR", "remotes"))
remotes::install_github("CBIIT/LDlinkR")
```

#### 1.2 Token

Make a one-time request for your personal access token from a web browser at [LDlink website](https://ldlink.nih.gov/?tab=apiaccess).

#### 1.3 Input files

This script utilizes one input file containing the rs identifiers of the markers preceded by the chromosome on which they are located:

1;rs1  
1;rs2  
1;rs3  
2;rs4  
3;rs5  
3;rs6  
4;rs7  
...

---

### 2. IBD and LR estimation 👨‍👩‍👧‍👦

This R script s designed to identify relationships between individuals based on their genotypic data.

#### 2.1 Install Required Packages

Make sure you have the required R packages installed. You can install them by running the following command in your R console:

```R
install.packages(c("pedtools", "devtools", "usethis"))
devtools::install_github("magnusdv/pedbuildr")
devtools::install_github("magnusdv/forrel")
```

#### 2.2 Input files

This script utilizes two input files:  
1.👫 A file containing **genetic information** of individuals, following the schema below:

| id  | fid | mid | sex | rs1 | rs2 | rs3 | rs4 | rs5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1   | 0   | 0   | 1   | 1/2 | 1/2 | 1/1 | 1/2 | 1/2 |
| 2   | 0   | 0   | 2   | 1/1 | 1/2 | 1/2 | 2/2 | 1/2 |
| 3   | 1   | 2   | 1   | 1/1 | 1/2 | 1/2 | 2/2 | 1/1 |

-   `id` : Individual ID 👤
-   `fid`: Father's ID 👨
-   `mid`: Mother's ID 👩
-   `sex`: Sex (1 for male, 2 for female) ⚤
-   `<rs1>`, `<rs2>`, `<rs3>`, etc.: Genotypes for each marker (e.g., 1/2, 1/1, 2/2) 🧬

    2.🧬 A file containing **information about markers** (RS IDs and allele frequencies), following the schema below:

rs1;  
1;0.5  
2;0.5  
rs2;  
1;0.75  
2;0.25  
rs3;  
1;0.33  
2;0.33  
3;0.33  
...

---

## Python 🐍 scripts

### Install Python

If you haven't already, you'll need to install a suitable Python3 distribution. You can download it from the [official Python website](https://www.python.org).

### 1. Install requirements

#### 1.1 Create a virtual environment

In a terminal window move to the desired folder and type:

`python -m venv .venv`

This will create a folder called `.venv` at the location where you executed the command.

#### 1.2 Activate the virtual environment

Execute the corresponding terminal command for your OS:

##### Linux/MacOS

`source .venv/bin/activate`

##### Windows

`.\venv\Scripts\activate.ps1`

#### 1.3 Install the requirements

Copy the `requirements.txt` file present in the Python folder to the working directory.

Run in the same terminal:

`pip install -r requirements.txt`

### 2. Execute the notebook code

Move the `kinship.ipynb` and `kinship_data.xlsx` files in the working directory.

Open it a in a Notebook compatible editor (e.g., VSCode, Spyder, etc.) and execute it to reproduce the experiments conducted.

Code is commented to facilitate its navigation.
