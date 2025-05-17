#' Run SuSiE Fine-Mapping Across Chromosomes Using MultiFun Prior and Summary Statistics
#'
#' This function runs SuSiE fine-mapping on one or all chromosomes using prior information from MultiFun,
#' GWAS summary statistics, LD matrices, and phenotype data. It wraps around \code{process_chr_prior}
#' for per-chromosome processing.
#'
#' @param Pheno_data A data frame containing phenotype data with at least two columns: ID and the phenotype variable.
#' @param Pheno_var Character string specifying the column name in \code{Pheno_data} that contains the phenotype.
#' @param Prior_path Base path to directory containing prior files (e.g., MultiFun probability files).
#' @param Prior_var Column name of the prior weight to use from prior files.
#' @param Summary_path Base path to directory containing GWAS summary statistic files.
#' @param LDMatrix_path Base path to directory containing LD matrix files.
#' @param Python_path Optional character string indicating the Python executable path needed to load sparse LD matrices.
#' @param Chr Character or numeric indicating which chromosome(s) to process. Use "All" to run all 22 autosomes.
#' @param Region Optional specification for a specific region to process; passed directly to \code{process_chr_prior}.
#' @param Output_path Path to output directory where results will be saved.
#'
#' @return A list of results returned by \code{process_chr_prior} for each processed chromosome.
#'
#' @export
#' @examples
#' \dontrun{
#' SuSiE_prior_summary(
#'   Pheno_data = my_pheno_df,
#'   Pheno_var = "TC",
#'   Prior_path = "/data/", #contain MultiFun folder
#'   Prior_var = "PriorCausalProbability_MKL",
#'   Summary_path = "/data/GWAS",
#'   LDMatrix_path = "data/UKBB_LD",
#'   Python_path = "/usr/bin/python3",
#'   Chr = "All",
#'   Output_path = "/results" #generate MultiFunSuSiE folder
#' )
#' }
SuSiE_prior_summary<-function( Pheno_data,
                               Pheno_var,
                               Prior_path,
                               Prior_var,
                               Summary_path,
                               LDMatrix_path,
                               Python_path=NULL,
                               Chr=22,  # Can also be "All"
                               Region="All",
                               Output_path){
  
  # Pheno_data=Pheno
  # Pheno_var=Phenotype
  # Prior_path=Prior_path
  # Prior_var=Prior_var
  # Summary_path=Summary_path
  # LDMatrix_path=LDMatrix_path
  # Python_path=NULL
  # Chr=22  #Chr="All"
  # Region="All"
  # Output_path=Output_path
  
  
  # Format Pheno
  Pheno <- base::subset(Pheno_data, select = c("ID", Pheno_var))
  base::colnames(Pheno)<-c("FID","Pheno")
  
  # Remove missing phenotypes
  if (any(is.na(Pheno$Pheno))) {
    Pheno <- Pheno[!is.na(Pheno$Pheno), , drop = FALSE]
  }
  
  # # Set Python path globally if provided
  # if (!is.null(Python_path)) {
  #   reticulate::use_python(Python_path, required = TRUE)
  # }
  
  # Helper function to create folder if not exists
  ensure_output_dir <- function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
    return(path)
  }
  Output_path <- ensure_output_dir(Output_path)
  
  # Run MultiFun+SuSiE for All or Specific Chromosome
  # Process chromosomes
  SuSiE_Result <- list()
  
  if (identical(Chr, "All")) {
    chr_list <- as.character(1:22)
    message("Running SuSiE for all 22 chromosomes...")
    
    SuSiE_Result <- lapply(chr_list, function(chr_list_tmp) {
      process_chr_prior(Pheno=Pheno,
                        Prior_path=Prior_path,
                        Prior_var=Prior_var,
                        Summary_path=Summary_path,
                        LDMatrix_path=LDMatrix_path,
                        Python_path=Python_path,
                        Chr_select=chr_list_tmp,
                        Region=Region,
                        Output_path=Output_path)
    })
    names(SuSiE_Result) <- paste0("CHR",chr_list)
  }
  
  if (!identical(Chr, "All")) {
    message(paste("Running SuSiE for chromosome:", Chr))
    
    SuSiE_Result[[paste0("CHR", Chr)]] <- process_chr_prior(Pheno=Pheno,
                                                            Prior_path=Prior_path,
                                                            Prior_var=Prior_var,
                                                            Summary_path=Summary_path,
                                                            LDMatrix_path=LDMatrix_path,
                                                            Python_path=Python_path,
                                                            Chr_select=Chr,
                                                            Region=Region,
                                                            Output_path=Output_path)
                                     
  }
  
  message("SuSiE analysis completed.")
  return(SuSiE_Result)
}


#' Process All or Specific Regions on a Single Chromosome Using SuSiE with MultiFun Prior
#'
#' Processes all defined regions on a chromosome in parallel using multiple cores.
#' It wraps around \code{process_region} for per-region processing, supports both Linux and Windows,
#' and handles failed runs via retry mechanism. Results are merged, deduplicated by position (POS),
#' and saved to CSV.
#'
#' @param Pheno A data frame containing phenotype information with columns "FID" and "Pheno".
#' @param Prior_path Base path to directory containing prior files (e.g., MultiFun probability files).
#' @param Prior_var Column name of the prior weight to use from prior files.
#' @param Summary_path Base path to directory containing GWAS summary statistic files.
#' @param LDMatrix_path Base path to directory containing LD matrix files.
#' @param Python_path Optional character string indicating the Python executable path needed to load sparse LD matrices.
#' @param Chr_select Character string specifying the chromosome number to process.
#' @param Region Character string specifying which region(s) to process; use "All" to run all regions.
#' @param Output_path Path to output directory where results will be saved.
#'
#' @return A data.table containing fine-mapping results for the processed chromosome.
#'         Returns NULL if merging fails.
#'
#' @export
#' @examples
#' \dontrun{
#' process_chr_prior(
#'   Pheno = my_pheno_df,
#'   Prior_path = "/data/", #contain MultiFun folder
#'   Prior_var = "PriorCausalProbability_MKL",
#'   Summary_path = "/data/GWAS",
#'   LDMatrix_path = "data/UKBB_LD",
#'   Python_path = "/usr/bin/python3",
#'   Chr_select = "1",
#'   Region = "All",
#'   Output_path = "/results" #generate MultiFunSuSiE folder
#' )
#' }
process_chr_prior <- function(Pheno,
                              Prior_path,
                              Prior_var,
                              Summary_path,
                              LDMatrix_path,
                              Python_path,
                              Chr_select,
                              Region,
                              Output_path) {
  
  # Load Region file
  Region_file_name <- paste0("Region_chr", Chr_select)
  Region_all <- get(Region_file_name, envir = parent.frame())
  
  # Determine Regions
  if (identical(Region, "All")) {
    Region_final<-seq_len(nrow(Region_all))
  }
  
  if (!identical(Region, "All")){
    #Region<-"15000001-54000007"
    Region_start_tmp<-as.numeric(strsplit(Region,'-')[[1]][1])
    Region_end_tmp<-as.numeric(strsplit(Region,'-')[[1]][2])
    if(Region_start_tmp<Region_all$start[1]){
      Region_start_row<-1
    }else{
      Region_start_row<-max(which(Region_all$start<=Region_start_tmp))
    }
    if(Region_end_tmp>Region_all$start[nrow(Region_all)]){
      Region_end_row<-nrow(Region_all)
    }else{
      Region_end_row<-min(which(Region_all$end>=Region_end_tmp))
    }
    
    Region_final<-c(Region_start_row:Region_end_row)
  }
  
  # Determine whether it is linux or windows
  os <- Sys.info()['sysname']
  
  if (os == "Linux") {
    cat("Running on Linux system.\n")
    
    # Determine number of cores
    num_cores <- if (!is.null(Sys.getenv("PBS_NUM_PPN"))) {
      as.integer(Sys.getenv("PBS_NUM_PPN"))
    } else {
      parallel::detectCores() - 10
    }
    num_cores <- max(5, num_cores) # Use at least 5 core
    
    # # # Run MultiFun+SuSiE for all selected Regions in parallel
    # Parallel processing with mclapply: It can only be applied to Linux systems and is not suitable for Windows systems
    SuSiE_CHR_All <- parallel::mclapply(Region_final, function(Region_select)
      process_region(Pheno=Pheno,
                     Prior_path=Prior_path,
                     Prior_var=Prior_var,
                     Summary_data=Summary_path,
                     LDMatrix_path=LDMatrix_path,
                     Python_path=Python_path,
                     Chr_select=Chr_select,
                     Region_select=Region_select),
      mc.cores = num_cores)
    
    # Retry loop for failed elements
    SuSiE_CHR_Valid <- SuSiE_CHR_All
    repeat {
      invalid_indices <- which(!sapply(SuSiE_CHR_Valid, inherits, what = "data.frame"))
      if (length(invalid_indices) == 0) {
        cat("All elements are now valid.\n")
        break  
      } else {
        cat("Retrying invalid regions:", paste(invalid_indices, collapse = ", "), "\n")
        reprocessed <- parallel::mclapply(Region_final, function(Region_select)
          process_region(Chr_select, Region_select, Region, Pheno, Prior_var, folder),
          mc.cores = num_cores)
        
        SuSiE_CHR_Valid[invalid_indices] <- reprocessed
      }
    }
    
  }else if (os == "Windows") {
    cat("Running on Windows system.\n")
    
    SuSiE_CHR_Valid <- lapply(Region_final, function(Region_select) {
      process_region(Pheno=Pheno,
                     Prior_path=Prior_path,
                     Prior_var=Prior_var,
                     Summary_data=Summary_path,
                     LDMatrix_path=LDMatrix_path,
                     Python_path=Python_path,
                     Chr_select=Chr_select,
                     Region_select=Region_select)
    })
    names(SuSiE_CHR_Valid) <- Region_final
  }
  
  # Merge results
  SuSiE_CHR_Result <- tryCatch({
    data.table::rbindlist(SuSiE_CHR_Valid, fill = TRUE)
  }, error = function(e) {
    warning("Merging results failed: ", e$message)
    return(NULL)
  })
  
  if (!is.null(SuSiE_CHR_Result)) {
    data.table::setkey(SuSiE_CHR_Result, POS)
    SuSiE_CHR_Result <- SuSiE_CHR_Result[, .SD[which.max(PIP)], by = POS]
  }
  
  # Save output
  out_folder <- file.path(Output_path, "MultiFunSuSiE")
  if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)
  
  filename <- paste0("MultiFunSuSiE_S_SNP_", gsub("PriorCausalProbability_", "", Prior_var), "_chr", Chr_select, ".csv")
  data.table::fwrite(SuSiE_CHR_Result, file.path(out_folder, filename))
  
  cat("Completed processing chromosome ", Chr_select, "\n")
  return(SuSiE_CHR_Result)
}


#' Process a Single Genomic Region Using SuSiE with MultiFun Prior
#'
#' This function processes a single genomic region by:
#' - Loading region and BIM files,
#' - Matching GWAS summary statistics with LD matrix SNPs,
#' - Applying MultiFun prior if available,
#' - Running SuSiE fine-mapping via \code{SuSiE_Summary},
#' - Returning posterior inclusion probabilities (PIP) for SNPs in the region.
#'
#' @param Pheno A data frame containing phenotype information with columns "FID" and "Pheno".
#' @param Prior_path Base path to directory containing prior files (e.g., MultiFun probability files).
#' @param Prior_var Column name of the prior weight to use from prior files.
#' @param Summary_data Base path to directory containing GWAS summary statistic files.
#' @param LDMatrix_path Base path to directory containing LD matrix files.
#' @param Python_path Optional character string indicating the Python executable path needed to load sparse LD matrices.
#' @param Chr_select Character string specifying the chromosome number to process.
#' @param Region_select Integer index indicating which region on the chromosome to process.
#'
#' @return A data.frame with two columns: 'POS' (position) and 'PIP' (posterior inclusion probability),
#'         or an empty data.frame if no valid data is found.
#'
#' @keywords internal
#' @importFrom data.table fread
#' @importFrom stats var
#' @importFrom utils modifyList
#' @importFrom susieR susie_rss
#' @importFrom reticulate use_python py_install import
#' @noRd
process_region <- function(Pheno,
                           Prior_path,
                           Prior_var,
                           Summary_data,
                           LDMatrix_path,
                           Python_path,
                           Chr_select,
                           Region_select) {
  
  # Region
  # Load Region file
  Region_file_name <- paste0("Region_chr", Chr_select)
  Region_all <- get(Region_file_name, envir = parent.frame())
  
  Region_start <- Region_all$start[Region_select]
  Region_end <- Region_all$end[Region_select]
  
  # Load BIM file
  bim_file_name <- paste0("ukb_bim_chr", Chr_select)
  #data(list = bim_file_name, package = "MultiFunPRS", envir = parent.frame())
  BIM <- get(bim_file_name, envir = parent.frame())
  
  # Load summary statistics
  gwas_file <- list.files(path = Summary_path, pattern = as.character(Chr_select), full.names = FALSE)
  if(length(gwas_file)!=1){
    stop(paste0("Not found Chr",Chr_select,"'s Summary statistics file or Found multiple file."))
  }
  gwas_file<-file.path(Summary_path, gwas_file)
  GWAS_summary <- data.table::fread(gwas_file, header = TRUE)
  colnames(GWAS_summary)[which(colnames(GWAS_summary) %in% "BP")]<-"POS"
  GWAS_summary <- subset(GWAS_summary, Region_start <= POS & POS <= Region_end)
  GWAS_summary <- subset(GWAS_summary, POS %in% BIM$POS)

  if (nrow(GWAS_summary) == 0) {
    return(data.frame(POS = numeric(), PIP = numeric()))
  }
  
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
  
  LD_file <- file.path(LDMatrix_path, paste0("chr", Chr_select) ,paste0("chr", Chr_select, "_", Region_start, "_", Region_end, ".npz"))
  LD_arr <- tryCatch({ load_sparse_LD_matrix(LD_file,Python_path) }, error = function(e) NULL)
  
  if (is.null(LD_arr)) {
    warning("Could not load LD matrix for Chr",Chr_select," region ", Region_start, "-", Region_end)
    return(data.frame(POS = GWAS_summary$POS, PIP = 0))
  }
  
  # Match LD SNPs with GWAS SNPs
  snp_file <- file.path(LDMatrix_path, paste0("chr", Chr_select), paste0("chr", Chr_select, "_", Region_start, "_", Region_end, ".gz"))
  LD_snp <- data.table::fread(snp_file, header = TRUE)
  colnames(LD_snp) <- c("rsID", "CHR", "POS", "A1", "A2")
  tmp <- GWAS_summary[, c("POS", "A1", "A2")]
  tmp$index <- 1
  LD_snp <- merge(LD_snp, tmp, by = c("POS", "A1", "A2"), all.x = TRUE)
  valid_idx <- which(!is.na(LD_snp$index))
  LD_snp <- LD_snp[valid_idx, ]
  LD_arr <- LD_arr[valid_idx, valid_idx]
  
  GWAS_summary <- GWAS_summary[match(LD_snp$POS, GWAS_summary$POS), ]
  
  # Load Prior
  prior_file <- file.path(Prior_path, "MultiFun", paste0("MultiFun_Probability_CHR", Chr_select, ".txt"))
  if (file.exists(prior_file)) {
    Prior <- data.table::fread(prior_file, header = TRUE)
    Prior <- as.data.frame(Prior)
    colnames(Prior)[which(colnames(Prior) %in% "BP")]<-"POS"
    Prior <- subset(Prior, select = c("POS", Prior_var))
    Prior <- Prior[match(LD_snp$POS, Prior$POS), Prior_var]
    Prior[is.na(Prior)] <- 0
    has_prior <- !all(Prior == 0)
  } else {
    has_prior <- FALSE
    warning("Missing MultiFun folder or 'MultiFun_Probability.txt' for Chr",Chr_select," in ",Prior_path)
  }
  
  # Error-handled SuSiE call
  SuSiE_wrapper <- function(...) {
    tryCatch({
      SuSiE_Summary(...)
    }, error = function(e) {
      if (grepl("missing value where TRUE/FALSE needed", e$message)) {
        warning("Caught specific error in region ", Region_select, ": ", e$message)
        return(data.frame(POS = LD_snp$POS, PIP = 0))
      } else {
        stop(e)
      }
    })
  }
  
  if (has_prior && nrow(GWAS_summary) > 1) {
    result <- SuSiE_wrapper(GWAS_summary, LD_arr, Pheno, Prior)
  } else {
    result <- data.frame(POS = LD_snp$POS, PIP = 0)
  }
  
  # Clean up memory
  rm(BIM, GWAS_summary, LD_arr, LD_snp, Prior, tmp, valid_idx)
  gc()
  
  return(result)
}



#' Run SuSiE with Summary Statistics and LD Matrix Using Prior Information
#'
#' Perform SuSiE fine-mapping using GWAS summary statistics and precomputed LD matrix.
#'
#' @param Summary A data frame with at least the following columns: beta, se, POS.
#' @param R_ref A reference LD matrix (p x p).
#' @param Pheno A data frame with a 'Pheno' column indicating the trait values.
#' @param Prior A vector of prior weights for each SNP.
#' @return A data frame with two columns: 'POS' and 'PIP' (Posterior Inclusion Probability).
#' @export
#' @examples
#' \dontrun{
#' SuSiE_Summary(Summary_data, LD_matrix, Pheno_data, Prior_weights)
#' }
SuSiE_Summary <- function(Summary, R_ref, Pheno, Prior) {
  
  N <- nrow(Pheno)
  fit_rss <- susieR::susie_rss(bhat = Summary$beta, shat = Summary$se,
                               n = N, R = R_ref,
                               var_y = var(Pheno$Pheno),
                               L = 10,
                               residual_variance = 1,
                               estimate_prior_variance = FALSE,
                               verbose = TRUE,
                               coverage = 0.9,
                               prior_weights = as.matrix(Prior))
  
  SuSiE_Summary_Result <- as.data.frame(fit_rss$pip)
  SuSiE_Summary_Result$POS <- Summary$POS
  colnames(SuSiE_Summary_Result) <- c("PIP", "POS")
  SuSiE_Summary_Result <- SuSiE_Summary_Result[, c("POS", "PIP")]
  
  return(SuSiE_Summary_Result)
}

#' Run SuSiE with Genotype Data Using Prior Information
#'
#' Perform SuSiE fine-mapping using genotype data and prior weights.
#'
#' @param Geno A matrix of genotype data (n x p), where n is sample size and p is number of SNPs.
#' @param Pheno A data frame containing phenotype data, must have column 'Pheno'.
#' @param Prior A vector or matrix of prior weights for each SNP.
#' @return A data frame with two columns: 'SNP' and 'PIP' (Posterior Inclusion Probability).
#' @export
#' @examples
#' \dontrun{
#' SuSiE_Geno(Geno_bed, Pheno_data, Prior_weights)
#' }
SuSiE_Geno <- function(Geno, Pheno, Prior) {
  
  fit <- susieR::susie(as.matrix(Geno), Pheno$Pheno,
                       L = 10,
                       estimate_residual_variance = FALSE,
                       estimate_prior_variance = FALSE,
                       verbose = TRUE,
                       coverage = 0.9,
                       prior_weights = as.matrix(Prior))
  
  SuSiE_Geno_Result <- as.data.frame(fit$pip)
  SuSiE_Geno_Result$SNP <- rownames(SuSiE_Geno_Result)
  colnames(SuSiE_Geno_Result) <- c("PIP", "SNP")
  SuSiE_Geno_Result <- SuSiE_Geno_Result[, c("SNP", "PIP")]
  
  return(SuSiE_Geno_Result)
}
