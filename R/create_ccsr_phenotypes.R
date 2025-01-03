#' Map source codes to CCSR categories
#'
#' This function maps the source procedure codes to CCSR categories and reshapes the data frame wider in preparation for the ProcWAS. The input data frame must have 4 columns with "id" representing the patient id, "vocabulary" representing one of "cpt_hcpcs" or "icd10", "code" representing the procedure code, and "date" representing the date the procedure was performed.
#'
#' @param source_codes_df Input data frame with source codes
#' @param mapping_file Mapping file to be used to map source icd10, cpt codes to ccsr categories. Default is the icd10_cpt4_source_to_ccsr mapping file from the package.
#' @return Data frame with CCSR categories mapped and reshaped wider
#' @export

create_ccsr_phenotypes <- function(source_codes_df, mapping_file = ProcWAS::icd10_cpt4_source_to_ccsr) {
    # generate error messages for incompatible inputs
    if(sum(names(source_codes_df) %in% c("id", "vocabulary", "code", "date"))!= 4) {
      stop("Must supply a data frame with 4 columns: 'id', 'vocabulary', 'code', 'date'")
    }
    if(!class(source_codes_df[["code"]]) %in% c("character","factor")) {
        stop("Please ensure character or factor representation for codes.")
    }
    if(!all(source_codes_df$vocabulary %in% c("cpt_hcpcs", "icd10"))) {
  stop("The 'vocabulary' column must contain one of 'cpt_hcpcs' or 'icd10'.")
    }
    
    # print the number of rows and unique patients in the input data
    message(glue::glue("{nrow(source_codes_df)} rows, {source_codes_df %>% distinct(id)  %>% nrow()} participant ids in input data"))

   # for rows with mappable procedures, pull in ccsr information
   ## inner join with ccsr map
    all_mapped_codes_df <- source_codes_df   %>% 
        inner_join(mapping_file, by = c("vocabulary", "code"))   %>% 
        select(id, prccsr, prccsr_description, date)
   
    ## keep only rows where unique ProcWAS codes were applied on distinct dates for a given pt
   mapped_codes_on_distinct_dates_df <- all_mapped_codes_df  %>% 
        group_by(id, prccsr)   %>%
        summarise(distinct_date_count = length(unique(date)))   %>%
        ungroup()   %>%
        mutate(is_any_proc = case_when(
            distinct_date_count >= 1 ~ TRUE,
            distinct_date_count < 1 ~ FALSE)
            )   %>%
        select(-distinct_date_count)  
    
    ## Reshape mapped codes from longer to wider df
    mapped_codes_wider_df  <-  pivot_wider(data = mapped_codes_on_distinct_dates_df,
                        names_from = prccsr,
                        values_from = is_any_proc,
                        values_fill = FALSE)
    
    message(glue::glue("{nrow(mapped_codes_wider_df)} participant ids had procedures mapped to CCSR codes on distinct days"))
    
    # Find ids without a mapped code
    unmapped_ids_df <- source_codes_df  %>% 
    anti_join(mapped_codes_wider_df, by = "id")   %>% 
    distinct(id)
    

    if (nrow(unmapped_ids_df) > 0) {
        message(glue::glue("There were {nrow(unmapped_ids_df)} ids present in the input data frame without a mapped code. These need to be added in the covariates data to be included in the ProcWAS."))
    } else {
        message("All ids had a procedure that were mapped to at least one CCSR code.")
    }

    return(mapped_codes_wider_df)
}