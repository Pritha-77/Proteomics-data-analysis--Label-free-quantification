# Proteomics-data-analysis--Label-free-quantification
# MaxQuant Proteomics Analysis Pipeline

This repository provides a clean folder structure and analysis workflow for processing proteomics data using **MaxQuant** and downstream statistical analysis in **R** (LFQ, PCA, DE analysis, volcano plots, etc.).

---

## 📁 Folder Structure

~/MaxQuant/

├── input/

│   ├── raw/        # all .raw files

│   ├── fasta/      # FASTA protein database files

│   └── mqpar/      # MQPar.xml parameter files

├── bin/            # MaxQuant executables / scripts

└── dotnet8/        # .NET 8 runtime required by MaxQuant


## Description of Folders

- **input/raw/**: Place all your raw MS files here. MaxQuant will process these files.  
- **input/fasta/**: Contains the FASTA database(s) used for protein identification.  
- **input/mqpar/**: Contains MQpar.xml parameter files specifying experimental design, enzyme settings, and quantification parameters.  
- **bin/**: Executable files or scripts to run MaxQuant from command line (optional).  
- **dotnet8/**: Required .NET runtime for MaxQuant.  

---

## Usage

1. Prepare your raw files, FASTA database, and MQpar.xml file in the corresponding folders.  
2. Run MaxQuant using the MQpar.xml file from the `input/mqpar/` folder.  
3. After MaxQuant finishes, use the R scripts in this repository to:
   - Filter and process LFQ intensities
   - Perform log2 transformation and missing value imputation
   - Generate PCA plots, correlation heatmaps
   - Conduct differential protein expression analysis using `limma`
   - Create volcano plots and save results 


## 📚 Dependencies

### Software
- **MaxQuant** (any recent version)
- **.NET 8** runtime (included in `dotnet8/`)

### R Packages
`tidyverse`, `limma`, `pheatmap`, `RColorBrewer`,  
`readxl`, `imputeLCMD`, `ggrepel`, `readr`














