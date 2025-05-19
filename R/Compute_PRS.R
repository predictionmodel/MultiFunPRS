#' Compute Polygenic Risk Score (PRS) Based on Potential Causal SNPs and Their Weights
#'
#' This function calculates the Polygenic Risk Score (PRS) using a set of potential causal SNPs 
#' and their corresponding effect sizes (beta values). Missing genotypes are imputed as 0. 
#' For each SNP, if the beta coefficient is positive, the score is calculated by multiplying 
#' the genotype dosage with beta; if negative, it uses the complementary allele dosage.
#'
#' @param snp_list Character vector of SNP IDs to be used in PRS calculation.
#' @param genotype A data frame containing genotype data with rows as individuals and columns as SNPs. 
#'                 Must include an "ID" column for individual identifiers.
#' @param beta_values A data frame containing SNP effect sizes with at least two columns: 
#'                    "rsID" (SNP ID) and "beta" (effect size).
#'
#' @return A list containing one element:
#'         - `PRS`: A data frame with two columns: `ID` and `PRS`, representing individual IDs and their computed PRS scores.
#'
#' @examples
#' # Example usage:
#' prs_result <- compute_prs(snp_list = c("rs123", "rs456"), 
#'                           genotype = Geno,
#'                           beta_values = Beta)
#'
#' @export
Compute_PRS <- function(snp_list, genotype, beta_values) {
  
  # Impute missing genotypes as 0
  genotype[is.na(genotype)] <- 0
  
  # Initialize PRS column
  genotype$PRS <- 0
  
  # Loop through each SNP in the list
  for (i in seq_along(snp_list)) {
    var <- snp_list[i]
    beta <- beta_values[beta_values$rsID == var, "beta"]
    
    # Only proceed if beta exists and is valid
    if (length(beta) > 0 && !is.na(beta)) {
      if (beta >= 0) {
        genotype$PRS <- genotype$PRS + genotype[[var]] * beta
      } else {
        genotype$PRS <- genotype$PRS + (2 - genotype[[var]]) * abs(beta)
      }
    }
  }
  
  # Return only ID and PRS
  result <- list(PRS = genotype[, c("ID", "PRS")])
  return(result)
}