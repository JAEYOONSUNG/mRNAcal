# mRNAcal
mRNAcal is an R package designed to streamline the analysis of mRNA regions and RNA secondary structures. This package offers tools for organizing GenBank data, extracting FASTA sequences, generating mRNA regions, and processing RNA secondary structures using viennaRNA package.

## Features
- GenBank Data Organization: Easily manage and extract information from GenBank files.
- FASTA Extraction: Extract coding sequences (CDS) or other genomic regions in FASTA format.
- Custom Downstream Sequences: Define a specific gene of interest (GOI) sequence to append to upstream regions.
- mRNA Region Generation: Automatically generate mRNA regions with user-defined upstream lengths and optional GOI.
- RNA Secondary Structure Prediction: Integrate RNAfold to calculate dot-bracket structures and minimum free energy (ΔG).
- Pipeline Automation: Use run_mRNAcal() to execute the full pipeline in one step with customizable parameters.


## Installation

To install the mRNAcal package and its dependencies, please follow these steps:
Check if devtools is installed, and install it if necessary
```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
```
### Step 1: Install Bioconductor dependencies
First, install the necessary Bioconductor packages:

```r
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "ComplexHeatmap", "Peptides"))
```

### Step 2: Install CRAN dependencies

Next, install the required CRAN packages. Some packages may not be installed automatically, so you can install them manually:
 sometimes it can not installed automatically followed

# Install CRAN packages
```r
install.packages(c("qdap", "seqinr", "splitstackshape"))
#remotes::install_github("trinker/qdap")
```

### Step 1: Install Java Development Kit (JDK)

The `rJava` package requires Java to be installed. Please install the Java Development Kit (JDK) from the following sources:

- **Windows/Mac**: [Oracle JDK Downloads](https://www.oracle.com/java/technologies/javase-jdk11-downloads.html)
- **Linux**: Use your package manager to install OpenJDK. For example:
  ```bash 
  sudo apt install openjdk-11-jdk

After installation, ensure that the JAVA_HOME environment variable is set correctly.

    Windows:

        1.    Go to System Properties -> Environment Variables -> System Variables -> New
        2.    Variable name: JAVA_HOME
        3.    Variable value: C:\Program Files\Java\jdk-11.0.1 (your JDK installation path)

    Mac/Linux:
        Add the following to your .bash_profile or .zshrc:
        ```bash
        export JAVA_HOME=$(/usr/libexec/java_home)

### Step 2: Install rJava

Once Java is installed and configured, you can install the rJava package:
    Windows:
        ```bash
        Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk-11.0.1")
        ```
        ```r
        install.packages("rJava")
        ```
    Mac/Linux:
        ```r
        install.packages("rJava")
        ```

### Step 3: Install mRNAcal package
devtools::install_github("JAEYOONSUNG/mRNAcal")

## Usage
1. Organizing GenBank Data
```r
# Organize GenBank file into a structured table
Genbank_organizer("example.gbk")
```

2. Extract FASTA Sequences
```r
genbank_fna_extractor()
```
3. Generate mRNA Regions
```r
generate_mRNA_regions(upstream_values = seq(100, 500, 100), downstream = 100)
# Define upstream and downstream regions of coding sequences for mRNA.
# In cases where you have a specific gene of interest (GOI) sequence to append to upstream regions, you can use the custom_downstream_seq parameter.
```
| locus_tag | seq_u100_dGOI                                | seq_u200_dGOI                                | ... |
|-----------|---------------------------------------------|---------------------------------------------|-----|
| gene1     | `ATCGG...ATGAAAAAAAAAAAAAAAAAAAAAAAAAAAA`   | `GCTAG...ATGAAAAAAAAAAAAAAAAAAAAAAAAAAAA`   | ... |
| gene2     | `CGTTA...ATGAAAAAAAAAAAAAAAAAAAAAAAAAAAA`   | `TACGG...ATGAAAAAAAAAAAAAAAAAAAAAAAAAAAA`   | ... |

---
4. Process RNAfold Results
The resulting `genbank_table` will include columns for each combination of upstream length and downstream region, as well as the custom sequence. Example:

# Process RNA sequences using RNAfold
process_all_sequences(data = "genbank_table", rnafold_path = "/path/to/RNAfold")

## Quick start
5. Automate the Entire Pipeline
#### native sequences
```r
run_mRNAcal(
  genbank_file = "example.gbk",
  rnafold_path = "/path/to/RNAfold",
  upstream_values = seq(100, 500, 100),
  downstream = 100
)
```
#### specific sequence
```r
run_mRNAcal(
  genbank_file = "example.gbk",
  rnafold_path = "/path/to/RNAfold",
  upstream_values = seq(100, 500, 100),
  custom_downstream_seq = "ATGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
```

## Result

| **Sequence**   | **Structure**                | **Free Energy (ΔG)** |
|----------------|------------------------------|----------------------|
| seq_u100_d100  | ..((((((((((.(((((((.......)))))...))))))..  | -51.3 kcal/mol |
| seq_u200_d100  | ....(((.(((.((((((((...).)))))))))))).)))..  | -85.6 kcal/mol |
| seq_u300_d100  | ...((((...((((.((((....))))..))))..........  | -102.4 kcal/mol |
| seq_u400_d100  | ((((((........((((((...)............)))))).  | -138.4 kcal/mol |
| seq_u500_d100  | ....((..((((.(((((((...))).))).))))..))....  | -39.3 kcal/mol |


### To be updated
- plot mRNA fold
- plot heatmap generation using free energy values
- delete duplicated mRNA sequences in excel file
- ciruclar form auto detection algorithm
