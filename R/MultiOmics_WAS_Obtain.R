#' Full MultiFun Framework Function
#'
#' Input phenotype and run full MultiFun framework to obtain multi-omics summary statistics.
#'
#' @param Phenotype Character string specifying the trait of interest.
#' @param folder Path to the working directory where output will be saved.
#' @return Returns path to processed GWAS file or creates necessary files in specified folder.
#' @export
#' @examples
#' \dontrun{
#' MultiOmics_WAS_Obtain(Phenotype = "AFib", folder = "/public/home/test_blank/tyd/Application/Result/AFib")
#' }

MultiOmics_WAS_Obtain<-function(Phenotype,folder){
  
  #(1.1) GWAS summary statistics:####
  # Set up GWAS output folder
  GWAS_folder<-paste0(folder,"/GWAS")
  
  # Create GWAS folder if it does not exist
  if(dir.exists(GWAS_folder)==F){
    system(paste0("mkdir ",GWAS_folder)) #Linux
    
    if(dir.exists(GWAS_folder)==F){
      dir.create(GWAS_folder)  #Windows
    }
    
    # Load phenotype ID mapping
    #data(list = "UKBB_94traits_release1", package = "MultiFunPRS", envir = parent.frame())
    Phenotype_ID_file<-get("UKBB_94traits_release1", envir = parent.frame())
    Phenotype_ID<-Phenotype_ID_file[which(Phenotype_ID_file$trait==Phenotype),"Phenotype_Code"]
    
    # Get corresponding GWAS file name
    #data(list = "UKBB_GWAS_Imputed_v3", package = "MultiFunPRS", envir = parent.frame())
    Phenotype_file<-get("UKBB_GWAS_Imputed_v3", envir = parent.frame())
    Phenotype_file<-Phenotype_file[which(Phenotype_file$Phenotype.Code==Phenotype_ID&Phenotype_file$Sex=="both_sexes"),"File"]
    Phenotype_outfile<-gsub("tsv.bgz","txt",Phenotype_file)
    
    ### Step 1: Create a munged summary statistics file in friendly format ###
    # (1) Split first column to extract CHR, POS, A1, A2
    system(paste0("zcat E:/Example_MultiFunPRS/UKBB_GWAS_Imputed_v3/",Phenotype_file," |
                  awk 'NR==1 {print \"CHR\\tPOS\\tA1\\tA2\\t\" $0; next}
                  NR>1 {split($1, arr, \":\"); print arr[1], arr[2], arr[3], arr[4], $0}' OFS=\"\\t\" > ",GWAS_folder,"/",Phenotype_outfile))
    #system(paste0("head ",Phenotype_ID,"_raw.gwas.imputed_v3.both_sexes.varorder.txt")) #But there is no SNP IDs in this file
    
    # (2) Split into 22 chromosome files
    system(paste0("awk 'NR==1 {header=$0; next} {print > \"",GWAS_folder,"/CHR_\" $1 \".txt\"}
                        END {for (i=1; i<=22; i++) system(\"echo \\\"\" header \"\\\" | cat - ",GWAS_folder,"/CHR_\" i \".txt > ",GWAS_folder,"/tmp && mv ",GWAS_folder,"/tmp ",GWAS_folder,"/CHR_\" i \".txt\")}' ",GWAS_folder,"/",Phenotype_outfile))
    
    ### Step 2: Merge with UKB BIM files to get SNP IDs ###
    read_and_merge<-function(i) {
      bim_file_name <- paste0("ukb_bim_chr", i)
      data(list = bim_file_name, package = "MultiFunPRS", envir = parent.frame())
      bim_file<-get(bim_file_name, envir = parent.frame())
      bim_file<-bim_file[,-which(colnames(bim_file) %in% c("A1", "A2"))]
      
      GWAS_file <- fread(paste0(GWAS_folder, "/CHR_", i, ".txt"), header = TRUE)
      Merge_file <- merge(bim_file, GWAS_file, by = c("CHR", "POS"),)
      return(Merge_file)
    }
    
    results <- lapply(1:22, read_and_merge) 
    All_file <- rbindlist(results, fill = TRUE) 
    
    # Save final merged GWAS file
    data.table::fwrite(All_file, paste0(GWAS_folder, "/GWAS_ALL.txt"), row.names = FALSE,sep=" ")
  }
  
  ######################################
  #(1.2) TWAS summary statistics:####
  # Generate TWAS folder
  TWAS_folder<-paste0(folder,"/TWAS") 
  # Create TWAS folder if it does not exist
  if(dir.exists(TWAS_folder)==F){system(paste0("mkdir ",TWAS_folder))}
  
  # Obtain TWAS for Phenotype
  system(paste0("awk '$12 == \"",Phenotype,"\"' /public/home/test_blank/tyd/Applicataion/Data/TWAS/UKBB_94traits_release1.1/release1.1/UKBB_94traits_release1.bed/UKBB_94traits_release1.bed > ",TWAS_folder,"/TWAS.txt"))
  
  # Format TWAS
  TWAS<-fread(paste0(TWAS_folder,"/TWAS.txt"),header=F) 
  coln<-fread("/public/home/test_blank/tyd/Applicataion/Data/TWAS/UKBB_94traits_release1.1/release1.1/UKBB_94traits_release1.cols",header = F)
  colnames(TWAS)<-coln$V1
  TWAS<-TWAS[which(TWAS$method=="SUSIE"),] 
  TWAS$chromosome<-sub("chr","",TWAS$chromosome)
  colnames(TWAS)[match(c("chromosome","rsid","allele1","allele2"),colnames(TWAS))]<-c("CHR","rsID","A2","A1")
  TWAS$POS<-as.numeric(unlist(strsplit(as.character(TWAS$variant),split=":"))[seq(2,dim(TWAS)[1]*4,4)]) 
  
  # Save final TWAS file
  fwrite(TWAS,paste0(TWAS_folder,"/TWAS.txt"),row.names=F)
  
  ######################################
  #(1.3) PWAS summary statistics:####
  # Generate PWAS folder
  PWAS_folder<-paste0(folder,"/PWAS") 
  if(dir.exists(PWAS_folder)==F){system(paste0("mkdir ",PWAS_folder))}
  
  # Obtain Protein
  #Input X:Protein: From /mnt/ndisk1/Student/tyd/UKBClinicalData/UKB_Protein
  setwd("/public/home/test_blank/tyd/Application/Data/UKB_Protein")
  Protein <- fread("Instance_0_Plasma_Protein_qc.txt", header = TRUE)
  Protein <- as.data.frame(Protein)
  #head(Protein)[,1:3]
  Protein$ID <- Protein$participant.eid
  Protein_Var <- colnames(Protein)[-which(colnames(Protein) %in% c("olink_instance_0.eid", "participant.eid", "ID"))]
  #length(Protein_Var) #2921 Proteins
  
  # Obtain Phenotype Data
  setwd("/public/home/test_blank/tyd/Application/Data/Phenotype")
  Phenotype_data <- fread("1.5UKBB_94traits_qc.csv", header = TRUE)
  Phenotype_data <- as.data.frame(Phenotype_data)
  Phenotype_data$ID <- Phenotype_data$ID_2024_03 # Because covariates use ID_2024_03
  Phenotype_data <- Phenotype_data[, c("ID", "ID_2021_08", Phenotype)]
  
  # Obtain Covariates
  # Input Corvariate:From /mnt/ndisk1/Student/tyd/UKBClinicalData/2024_03
  setwd("/public/home/test_blank/tyd/Application/Data/UKB_Protein/Corvariate_protein")
  Corvariate_protein <- fread("1.5Corvariate_protein_qc.csv", header = TRUE)
  Corvariate_protein <- as.data.frame(Corvariate_protein)
  
  # Join datasets
  Data <- left_join(Protein, Phenotype_data, by = "ID")
  Data <- left_join(Data, Corvariate_protein, by = "ID")
  
  corvariate <- c("Age", "Sex", "BMI", "batch", "UKB_centre", "UKB_array_type", "UKB_PPP_subcohort", paste0("PC", 1:20), "Time_sample_measurement")
  
  linear_reg <- function(x) {
    formula <- as.formula(paste(Phenotype, "~ x +", paste(corvariate, collapse = "+")))
    summary <- summary(lm(formula, data = Data))
    return(coef(summary)[2, ])
  }
  
  logistic_reg <- function(x) {
    formula <- as.formula(paste(Phenotype, "~ x +", paste(corvariate, collapse = "+")))
    summary <- summary(glm(formula, data = Data, family = binomial))
    return(coef(summary)[2, ])
  }
  
  # Determine if Phenotype is binary and choose appropriate model
  clean_variable <- na.omit(Data[, Phenotype])
  is_binary <- length(unique(clean_variable)) == 2
  
  if(is_binary){
    Results <- as.data.frame(t(apply(Data[, Protein_Var], 2, logistic_reg)))
  } else {
    Results <- as.data.frame(t(apply(Data[, Protein_Var], 2, linear_reg)))
  }
  Results$Protein <- row.names(Results)
  fwrite(Results, paste0(PWAS_folder, "/Protein.csv")) # Generate Protein associated with Phenotype
  
  
  # Obtain PWAS for Phenotype
  # Obtain the UKBPPP.ProteinID and SNP corresponding to each Protein based on ST3 in Plasma proteomic associations with genetics and health in the UK Biobank
  Protein_ID <- read.csv("/public/home/test_blank/tyd/Application/Data/UKB_Protein/Protein_ID.csv", header = TRUE)
  Protein <- left_join(Results, Protein_ID, by = "Protein")
  Protein$Phenotype <- Phenotype
  colnames(Protein) <- c("Estimate", "Se", "t_value", "P_value", "Protein", "UKBPPP.ProteinID", "Olink.ID", "Assay.Target", "Phenotype")
  Protein <- Protein[which(Protein[, "P_value"] < (1.7*10^-5)), ] # Filter significant proteins
  Protein_SNP <- read.csv("/public/home/test_blank/tyd/Application/Data/UKB_Protein/Protein_SNP.csv", header = TRUE)
  PWAS <- left_join(Protein_SNP, Protein, by = "UKBPPP.ProteinID")
  PWAS <- PWAS[which(PWAS$Phenotype == Phenotype), ]
  colnames(PWAS)[1] <- "Variant_ID"
  PWAS <- PWAS[!duplicated(PWAS[, "Variant_ID"]), ]
  PWAS$CHR <- PWAS$CHROM
  
  # Save final PWAS file
  fwrite(PWAS, paste0(PWAS_folder, "/PWAS.txt"), row.names = FALSE)
}