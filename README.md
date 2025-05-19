# MultiFunPRS

# MultiFun and MultiFun-PRS
**MultiFun** (Multi-omics Functional Prediction)

**MultiFun-PRS** (Multi-omics Functional Polygenic risk score)


The **MultiFunPRS R package** contains the code of the methods MultiFun for causal SNPs prediction leveraging multi-moics functional annotations and summary statistics, MultiFun-PRS for polygeinc risk score with functionally-informed fine-mapping and MultiFun as prior. MultiFun and MultiFun-PRS are described in the Article "Leveraging multi-omics and functional annotations to identify causal SNPs and improve polygenic risk scores of complex traits within and across ancestries."

**MultiFun** estimates prior causal probabilities for SNPs,  which can then be used by fine-mapping methods like SuSiE. Unlike previous methods for functionally-informed fine-mapping, MultiFun can aggregate multi-omics layers data from across the entire genome and two large-scale functional annotation resources.

**MultiFun-PRS** exploits causal probabilities and fine-mapping to improve within and cross-ancestries polygenic risk scores, by predicting using causal effect estimates intead of tagging effect estimates.

# Manual
We provide a detailed manual of MultiFun and MultiFun-PRS in the Wiki page, including downloads of various resources. If you encounter any questions, please check the FAQ first.

# Installation



# Example

We provide an example to get up and running quickly.
Note: You need to download the "Example_MultiFunPRS" folder from Github in advance, which includes the summary and annotation data of the sample downloaded in advance. See the wiki page for the full data download.

```R
library(devtools)

# Set working directory to the root of your package
setwd("E:/MultiFunPRS")

# Build and install the package (optional --preclean for cache cleanup)
install(preclean = TRUE)

library(MultiFunPRS)

# View all datasets in the package
data(package = "MultiFunPRS")

# Run MultiFun ####
# 1. Input Multi-Omics Summary Statistics ####

# Approach 1: Use provided phenotype data
phenotype <- "TC"
folder <- "E:/Example_MultiFunPRS/TC"
chr_select <- 20

# 1) Obtain Multi-Omics summary statistics
# Note: This function can only be executed on Linux systems
# MultiOmics_WAS_Obtain(phenotype, Output_folder)

# 2) Format Multi-Omics summary statistics (e.g., GWAS/TWAS/PWAS):
MultiOmics_summary<-MultiOmics_WAS_format(Chr_select=Chr_select, # Use "All" to run all chromosomes
                                          folder=folder)

# Approach 2: Use own file###
## Users can also generate custom summary statistics from their own multi-omics data. Please ensure the following:
## A separate summary statistics file is required for each omics layer.
## Each file must contain the following columns: "CHR", "rsID", "CM", "POS", "A1", "A2", and a column representing statistical significance("Index"), such as P values or posterior inclusion probabilities (PIP) from fine-mapping. Alternatively, binary indicators (e.g., 0/1) can be used to denote significance based on other criteria.

# 2. Input Functional Annotations####
# Approach 1: Use Functional Annotations sources which we provided(FAVOR+BaselineLFmodel)###
# Need to download in advance, about 8.3G
Annotation_path <- paste0("E:/Example_MultiFunPRS/Annoatation/BaselineLF_FAVOR/Annot_BaselineLF_FAVOR_qc_", Chr_select, ".csv")
Annotation <- data.table::fread(Annotation_path, header = TRUE)
Annotation <- as.data.frame(Annotation)

# Approach 2: Use own file###

# 3. Run MultiFunR to obtain MultiFun Probability####
# For running on all chromosomes, Linux is recommended; Windows can also run but may take longer
MultiFun_Prior_Result<-MultiFun_Prior(Chr_select=Chr_select, # Use "All" to run all chromosomes
                                      folder,MultiOmics_summary,Annotation)

# Output: Create an output folder
MultiFun_folder <- file.path(folder, "/MultiFun")
if (!dir.exists(MultiFun_folder)) dir.create(MultiFun_folder)

# Save results
data.table::fwrite(MultiFun_Prior_Result$MultiFun_Weight, file.path(MultiFun_folder, "MultiFun_Weight.csv"))
MultiFun_Probability<-MultiFun_Prior_Result$MultiFun_Probability[[paste0("CHR",Chr_select)]]
data.table::fwrite(MultiFun_Probability, file.path(MultiFun_folder, paste0("MultiFun_Probability_CHR", Chr_select, ".txt")), row.names = FALSE)

################################################################

# Run MultiFun+SuSiE######
# 1. Input Phenotype####
Pheno_path<-file.path(folder, "Phenotype_TC.csv")
Pheno<-data.table::fread(Pheno_path,header=T)
Pheno<-as.data.frame(Pheno)
Pheno<-Pheno[,c("ID",Phenotype)]

# 2. Input Prior####
Prior_path <- folder
MultiFun_Prior<-c("PriorCausalProbability_MKL", "PriorCausalProbability_MAX")
Prior_var<-"PriorCausalProbability_MAX"

# 3. Input Summary statistics####
Summary_path<- file.path(folder, "GWAS")

#4. Input LD matrix####
LDMatrix_path<-"E:/Example_MultiFunPRS/UKBB_LD"
#Python_path<-"/public/software/miniconda3/bin/python3.9"

# Load sparse LD matrix
load_sparse_LD_matrix <- function(LD_file,Python_path) {
  if(!is.null(Python_path)){
    reticulate::use_python(Python_path, required = TRUE)
  }else{
    reticulate::use_python(reticulate::py_config()$python, required = TRUE)
  }
  cipy_available <- reticulate::py_module_available("scipy")
  if(!cipy_available){
    reticulate::py_install("scipy")
  }
  sparse <- reticulate::import("scipy.sparse")
  LD_arr <- as.array(sparse$load_npz(LD_file))
  return(LD_arr)
}

# 5. Define Output path ####
Output_path<-folder

# 6. Run functional-informed SuSiE to obtain PIP ####
# For running on all chromosomes, Linux is recommended for better performance
SuSiE_Pior_Result<-SuSiE_prior_summary(Pheno_data=Pheno,
                                       Pheno_var=Phenotype,
                                       Prior_path=Prior_path,
                                       Prior_var=Prior_var,
                                       Summary_path=Summary_path,
                                       LDMatrix_path=LDMatrix_path,
                                       Python_path=NULL,
                                       Chr=Chr_select, # Use "All" to run all chromosomes
                                       Region="All",
                                       Output_path=Output_path)
              



#################################################################

# Construct MultiFun-PRS
#1. Input Genotype of popential causal SNPs
library(plink2R)

# Suppose that we have screened causal SNPs based on PIP and obtained the genotypes of causal SNPs in the population of interest
TargetPopulation<-"White_British"
PRS_folder<-file.path(folder,"PRS") 
Geno_file <- file.path(PRS_folder,paste0("ukb_imp_chr", Chr_select,"_v3_qc_",TargetPopulation,"_Identified"))
Geno<-read_plink(Geno_file)
Geno<-as.data.frame(Geno$bed)
Identified_SNP<-colnames(Geno)
Geno$ID<-sub(":.*", "",row.names(Geno))


#2. Input Estimate of popential causal SNPs for weight
library(data.table)
library(dplyr)

# Load bim file
bim_file_name <- paste0("ukb_bim_chr", Chr_select)
BIM <- get(bim_file_name, envir = parent.frame())

Beta<-fread(paste0(folder,"/GWAS/CHR_",Chr_select,".txt"),header=T)
Beta<-as.data.frame(Beta)
Beta<-left_join(BIM,Beta,by=c("CHR","POS"))
Beta<-Beta[match(Identified_SNP,Beta$rsID),]

#3. Construct PRS
PRS_Result<-NULL
PRS_Result<-Compute_PRS(SNP=Identified_SNP,Geno=Geno,Beta=Beta)
names(PRS_Result)<-paste0("CHR",Chr_select)


```

