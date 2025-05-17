#' Format Multi-Omics WAS Data for MultiFun Framework
#'
#' This function formats GWAS, TWAS, and PWAS data for a specific chromosome to be used in the MultiFun framework.
#'
#' @param Chr_select Character or numeric indicating the chromosome number (e.g., "1", 1), or "All" to process all chromosomes.
#' @param GWAS_summary Data frame containing GWAS summary statistics with columns: CHR, POS, pval.
#' @param TWAS_summary Data frame containing TWAS summary statistics with columns: CHR, POS, method.
#' @param PWAS_summary Data frame containing PWAS summary statistics with columns: CHR, rsID, Phenotype.
#' @return A list of formatted data frames for each omic type per chromosome.
#' @export
#' @examples
#' \dontrun{
#' MultiOmics_WAS_format("All", folder)
#' MultiOmics_WAS_format("1", folder)
#' }

MultiOmics_WAS_format<-function(Chr_select,folder){

  if(Chr_select=="All"){

    MultiOmics_summary<-NULL
    for(All_chr_tmp in 1:22){
      # GWAS:
      GWAS_summary <- data.table::fread(paste0(folder, "/GWAS/CHR_", All_chr_tmp, ".txt"), header = TRUE)
      GWAS_summary <- as.data.frame(GWAS_summary)

      # TWAS
      TWAS_summary <- data.table::fread(paste0(folder, "/TWAS/TWAS.txt"), header = TRUE)
      TWAS_summary <- as.data.frame(TWAS_summary)

      # PWAS
      PWAS_summary <- data.table::fread(paste0(folder, "/PWAS/PWAS.txt"), header = TRUE)
      PWAS_summary <- as.data.frame(PWAS_summary)

      #Format MultiOmics_WAS by Function: MultiOmics_WAS_format
      MultiOmics_summary_tmp <- MultiOmics_WAS_format_chr(All_chr_tmp, GWAS_summary, TWAS_summary, PWAS_summary)
      rm(GWAS_summary, TWAS_summary, PWAS_summary)
      MultiOmics_summary[[All_chr_tmp]]<-MultiOmics_summary_tmp
      names(MultiOmics_summary)[All_chr_tmp]<-paste0("CHR",All_chr_tmp)
    }

  }else{

    MultiOmics_summary<-NULL

    # GWAS for chr:
    GWAS_summary <- data.table::fread(paste0(folder, "/GWAS/CHR_", Chr_select, ".txt"), header = TRUE)
    GWAS_summary <- as.data.frame(GWAS_summary)

    # TWAS
    TWAS_summary <- data.table::fread(paste0(folder, "/TWAS/TWAS.txt"), header = TRUE)
    TWAS_summary <- as.data.frame(TWAS_summary)

    # PWAS
    PWAS_summary <- data.table::fread(paste0(folder, "/PWAS/PWAS.txt"), header = TRUE)
    PWAS_summary <- as.data.frame(PWAS_summary)

    MultiOmics_summary[[1]] <- MultiOmics_WAS_format_chr(Chr_select, GWAS_summary, TWAS_summary, PWAS_summary)
    names(MultiOmics_summary)[1]<-paste0("CHR",Chr_select)
  }

  return(MultiOmics_summary)
}


#' Internal Function to Format Multi-Omics Data for a Single Chromosome
#'
#' Formats GWAS, TWAS, and PWAS data for a single chromosome and returns them in a list.
#'
#' @param Chr_select Character indicating the chromosome number (e.g., "1").
#' @param GWAS_summary Data frame containing GWAS summary statistics with columns: CHR, POS, pval.
#' @param TWAS_summary Data frame containing TWAS summary statistics with columns: CHR, POS, method.
#' @param PWAS_summary Data frame containing PWAS summary statistics with columns: CHR, rsID, Phenotype.
#' @return A list containing formatted GWAS, TWAS, and PWAS data frames for the specified chromosome.
#' @keywords internal
#' @importFrom dplyr left_join
#' @noRd
MultiOmics_WAS_format_chr<-function(Chr_select,GWAS_summary,TWAS_summary,PWAS_summary){

  # Load BIM data from package resources
  bim_file_name <- paste0("ukb_bim_chr", Chr_select)
  #data(list = bim_file_name, package = "MultiFunPRS", envir = parent.frame())
  BIM <- get(bim_file_name, envir = parent.frame())
  BIM <- BIM[, c("CHR", "rsID", "POS")]

  # GWAS processing
  GWAS_summary_Chr <- subset(GWAS_summary, CHR == Chr_select)
  GWAS_summary_Chr <- dplyr::left_join(BIM, GWAS_summary_Chr, by = c("CHR", "POS"))
  GWAS_summary_Chr$Index <- ifelse(GWAS_summary_Chr$pval < 5e-8, 1, 0)
  GWAS_summary_Chr[which(is.na(GWAS_summary_Chr$Index)),"Index"]<-0

  if(sum(GWAS_summary_Chr$Index)==0){
    GWAS_summary_Chr$Index <- ifelse(GWAS_summary_Chr$pval < 5e-6, 1, 0)
    GWAS_summary_Chr[which(is.na(GWAS_summary_Chr$Index)),"Index"]<-0
  }
  #table(GWAS_summary_Chr$Index)
  GWAS_Index_n <- sum(GWAS_summary_Chr$Index)
  GWAS_summary_Chr <- GWAS_summary_Chr[, c("CHR", "POS", "Index")]
  GWAS_summary_Chr <- GWAS_summary_Chr[!duplicated(GWAS_summary_Chr$POS), , drop = FALSE]

  # TWAS processing
  TWAS_summary_Chr <- subset(TWAS_summary, CHR == Chr_select)
  TWAS_summary_Chr <- dplyr::left_join(BIM, TWAS_summary_Chr, by = c("CHR", "POS"))
  TWAS_summary_Chr$Index <- ifelse(TWAS_summary_Chr$method == "SUSIE", 1, 0)
  TWAS_summary_Chr[which(is.na(TWAS_summary_Chr$Index)),"Index"]<-0
  TWAS_Index_n <- sum(TWAS_summary_Chr$Index)
  TWAS_summary_Chr <- TWAS_summary_Chr[, c("CHR", "POS", "Index")]
  TWAS_summary_Chr <- TWAS_summary_Chr[!duplicated(TWAS_summary_Chr$POS), , drop = FALSE]

  # PWAS processing
  PWAS_summary_Chr <- subset(PWAS_summary, CHR == Chr_select)
  PWAS_summary_Chr <- dplyr::left_join(BIM, PWAS_summary_Chr, by = c("CHR", "rsID"))
  PWAS_summary_Chr$Index <- ifelse(is.na(PWAS_summary_Chr$Phenotype), 0, 1)
  PWAS_Index_n <- sum(PWAS_summary_Chr$Index)
  PWAS_summary_Chr <- PWAS_summary_Chr[, c("CHR", "POS", "Index")]
  PWAS_summary_Chr <- PWAS_summary_Chr[!duplicated(PWAS_summary_Chr$POS), , drop = FALSE]

  # Combine results into list
  MultiOmics_summary <- list(
    GWAS = GWAS_summary_Chr,
    TWAS = TWAS_summary_Chr,
    PWAS = PWAS_summary_Chr
  )

  if (GWAS_Index_n == 0) MultiOmics_summary <- MultiOmics_summary[-which(names(MultiOmics_summary) == "GWAS")]
  if (TWAS_Index_n == 0) MultiOmics_summary <- MultiOmics_summary[-which(names(MultiOmics_summary) == "TWAS")]
  if (PWAS_Index_n == 0) MultiOmics_summary <- MultiOmics_summary[-which(names(MultiOmics_summary) == "PWAS")]

  return(MultiOmics_summary)
}

