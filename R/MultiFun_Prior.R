#' Process MultiOmics data for a single chromosome to calculate prior probabilities.
#'
#' This function processes one chromosome of multi-omics summary statistics and annotation data,
#' runs the MultiFun algorithm to compute prior causal SNP probabilities, and returns results.
#'
#' @param All_chr_tmp Chromosome number (integer).
#' @param MultiOmics_summary A list of multi-omics summary statistics by chromosome.
#' @param Annotation A data frame containing annotation features with at least columns: CHR, POS, and other features.
#' @return A list containing:
#' \item{MultiOmics_Weight}{Weights assigned to each omics layer.}
#' \item{MultiOmics_PriorProbability}{Data frame of predicted probabilities per omics and combined scores.}
#' @importFrom dplyr select
#' @export
Process_MultiFun_chr <- function(All_chr_tmp, MultiOmics_summary, Annotation) {

  # 1. Input Multi-omics summary statistics for this chromosome
  chr_key <- paste0("CHR", All_chr_tmp)
  MultiOmics_summary_chr <- MultiOmics_summary[[chr_key]]

  if (!is.null(MultiOmics_summary_chr)) {
    for(num in seq_along(MultiOmics_summary_chr)){
      colnames(MultiOmics_summary_chr[[num]])[colnames(MultiOmics_summary_chr[[num]]) == "BP"] <- "POS"
    }
  }

  # 2. Input Annotation data
  Annotation <- subset(Annotation, CHR == All_chr_tmp)
  colnames(Annotation)[colnames(Annotation) == "BP"] <- "POS"

  # 3. Run MultiFunR to obtain MultiFun Probability
  Result_list <- NULL
  if (!is.null(MultiOmics_summary_chr) && length(MultiOmics_summary_chr) > 0) {
    Result_list <- MultiFunR(MultiOmics_summary_chr, Annotation, "randomforest")
    message(paste0("MultiFun for CHR ", All_chr_tmp, " is completed."))
  } else {
    warning(paste0("MultiFun for CHR ", All_chr_tmp, " lacks MultiOmics-WAS."))
  }

  return(Result_list)
}

#' Main function to run MultiFun across selected chromosomes and generate prior probability scores.
#'
#' This function coordinates processing across multiple chromosomes, supports parallel execution on Linux,
#' and returns final weights and prior probabilities.
#'
#' @param Chr_select Either "All" to process all chromosomes or an integer indicating specific chromosome.
#' @param folder Path to helper functions (not used here but preserved for backward compatibility).
#' @param MultiOmics_summary A list of multi-omics summary statistics indexed by chromosome.
#' @param Annotation A data frame containing annotation features including CHR and POS.
#' @return A list containing:
#' \item{MultiFun_Weight}{A data frame of omics weights per chromosome.}
#' \item{MultiFun_Probability}{A named list of prior probability data frames per chromosome.}
#' @export
MultiFun_Prior <- function(Chr_select, folder, MultiOmics_summary, Annotation) {

  #1. Load helper functions
  # source('MultiFun_Function.R')

  # 2. Run MultiFun
  if (identical(Chr_select, "All")) {

    # Determine whether it is linux or windows
    os <- Sys.info()['sysname']
    MultiFun_Valid <- list()

    if (os == "Linux") {
      cat("Running on Linux system.\n")

      # Determine number of cores
      num_cores <- if (!is.null(Sys.getenv("PBS_NUM_PPN"))) {
        as.integer(Sys.getenv("PBS_NUM_PPN"))
      } else {
        parallel::detectCores() - 10
      }
      num_cores <- max(5, num_cores) # Use at least 5 core

      # # Process all chromosomes in parallel
      # Parallel processing with mclapply: It can only be applied to Linux systems and is not suitable for Windows systems
      MultiFun_All <- parallel::mclapply(1:22, function(All_chr_tmp)
        Process_MultiFun_chr(All_chr_tmp,MultiOmics_summary,Annotation),
        mc.cores = num_cores)

      # Validate and re-process invalid results
      # Due to the nuclear issue, the parallel operation keeps having problems, but no error is reported. It only leads to the inability to output MultiFun
      MultiFun_Valid <- MultiFun_All  # Initialize the MultiFun_Valid list

      # Use a while loop until all elements are valid
      repeat {
        invalid_indices <- which(!sapply(MultiFun_Valid, function(x) inherits(x, "list") || is.null(x))) # Find the element indexes that do not meet the conditions

        if (length(invalid_indices) == 0) {
          cat("All CHR are now valid.\n")
          break  # If there are no invalid elements, exit the loop
        } else {
          cat("Invalid CHR found at indices:", paste(invalid_indices, collapse = ", "), "\n")
          MultiFun_Reprocessed <- parallel::mclapply(1:22, function(All_chr_tmp)
            Process_MultiFun_chr(All_chr_tmp,MultiOmics_summary,Annotation),
            mc.cores = num_cores) # Reprocess the invalid elements

          MultiFun_Valid <- replace(MultiFun_Valid, invalid_indices, MultiFun_Reprocessed) #Put the reprocessed result into the new list
        }
      }
      names(MultiFun_Valid) <- paste0("CHR", 1:22)

    } else if (os == "Windows") {
      cat("Running on Windows system.\n")

      MultiFun_Valid <- lapply(X = 1:22, FUN = Process_MultiFun_chr,
                               MultiOmics_summary = MultiOmics_summary,
                               Annotation = Annotation)
      names(MultiFun_Valid) <- paste0("CHR", 1:22)

    } else {
      stop("Unsupported operating system: ", os)
    }

    # 3. Format and output results
    MultiFun_Weight <- NULL
    MultiFun_Probability <- list()

    for (chr_idx in 1:22) {
      chr_key <- paste0("CHR", chr_idx)
      result <- MultiFun_Valid[[chr_key]]

      if (!is.null(result) && !is.null(result$MultiOmics_Weight)) {
        weight_df <- as.data.frame(t(result$MultiOmics_Weight))
        prob_df <- result$MultiOmics_PriorProbability

        colnames(weight_df) <- colnames(prob_df)[grep("WAS", colnames(prob_df), value = TRUE)]
        weight_df$CHR <- chr_idx
        MultiFun_Weight <- dplyr::bind_rows(MultiFun_Weight, weight_df)

        MultiFun_Probability[[chr_key]] <- prob_df
      }
    }

    return(list(
      MultiFun_Weight = MultiFun_Weight,
      MultiFun_Probability = MultiFun_Probability
    ))

  } else {
    # Process only selected chromosome
    chr_num <- as.integer(Chr_select)
    result <- Process_MultiFun_chr(chr_num, MultiOmics_summary, Annotation)

    MultiFun_Weight <- NULL
    MultiFun_Probability <- list()

    if (!is.null(result)) {
      MultiFun_Weight <- as.data.frame(t(result$MultiOmics_Weight))
      MultiFun_Probability[[paste0("CHR", chr_num)]] <- result$MultiOmics_PriorProbability

      colnames(MultiFun_Weight) <- colnames(MultiFun_Probability)[grep("WAS", colnames(MultiFun_Probability))]
      MultiFun_Weight$CHR <- chr_num
    }

    return(list(
      MultiFun_Weight = MultiFun_Weight,
      MultiFun_Probability = MultiFun_Probability
    ))
  }
}
