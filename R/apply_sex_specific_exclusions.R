#' Apply sex specific exclusions
#'
#' This function applies sex-specific exclusions for male and female reproductive system procedures
#' FRS___ = female reproductive system
#' MRS___ = male reproductive system
#' PGN___ = pregnancy related procedures
#'
#' @param phewas_df Input data frame
#' @param name_of_sex_column Name of the sex column in the input data frame
#' @return phewas data frame with sex exclusions applied
#' @export

apply_sex_specific_exclusions <- function(phewas_df, name_of_sex_column){
    # Check if the sex column contains only valid values
  valid_sex_values <- c("M", "m", "F", "f", "Male", "male", "Female", "female", "Intersex", "intersex", "None", "none", "Other", "other")
  if (!all(phewas_df[[name_of_sex_column]] %in% valid_sex_values)) {
    stop("The sex column contains invalid values. Must be one of 'M', 'm', 'Male, 'male', 'F', 'f', 'Female', 'female', 'Intersex', 'intersex', 'None', 'none', 'Other', 'other'")
  }
  
  phewas_sex_exclusions_df <- phewas_df %>% 
    # if code is female specific and sex is male then convert to NA
    mutate(across(starts_with("PGN") | starts_with("FRS"), 
                  ~ if_else(!!sym(name_of_sex_column) %in% c("M", "m", "Male", "male"), NA, .)))  %>% 
    # if code is male specific and sex is female then convert to NA
    mutate(across(starts_with("MRS"), 
                  ~ if_else(!!sym(name_of_sex_column) %in% c("F", "f", "Female", "female"), NA, .)))  %>% 
  
  return(phewas_sex_exclusions_df)
}