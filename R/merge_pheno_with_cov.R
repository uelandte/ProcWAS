#' Merge phenotypes with covariates
#'
#' This function merges the phenotype information (output of create_ccsr_phenotypes) with covariate information.
#'
#' @param phenotypes_df Input phenotypes data frame
#' @param covariates_df Input covariates data frame
#' @param append_cov_without_pheno_as_controls Should individuals in the covariate data but not in the phenotypes data be appended as controls to the phenotypes data Default is TRUE
#' @return Merged phenotypes and covariates data frame
#' @export

merge_pheno_with_cov <- function(phenotypes_df, covariates_df, append_cov_without_pheno_as_controls = TRUE){

cov_but_no_pheno_info_df <- covariates_df  %>% 
    anti_join(phenotypes_df, by = "id")  %>% 
    distinct(id)
    
message(glue::glue("{nrow(cov_but_no_pheno_info_df)} rows in covariates data do not have a match in phenotypes data."))

    if (nrow(cov_but_no_pheno_info_df) == 0) {
        pheno_cov_df <- covariates_df  %>% 
            inner_join(phenotypes_df, by = "id")  %>% 
            mutate(across(where(is.logical), ~replace_na(., FALSE)))
    } else if (nrow(cov_but_no_pheno_info_df) > 0) {
        if(append_cov_without_pheno_as_controls == TRUE) {
            message("Appending these individuals as controls to the input data.")
            pheno_cov_df <- covariates_df  %>% 
                left_join(phenotypes_df, by = "id")  %>% 
                mutate(across(where(is.logical), ~replace_na(., FALSE))) # if not in pheno file, add FALSE for the proc codes 
        } else {
            message("These individuals will be excluded from the phenotypes data.")
            pheno_cov_df <- covariates_df  %>% 
                inner_join(phenotypes_df, by = "id")  %>% 
                mutate(across(where(is.logical), ~replace_na(., FALSE))) 
        }
    }
return(pheno_cov_df)
}