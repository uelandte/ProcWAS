# Set seed for reproducibility
set.seed(123)

# Create source codes data frame
source_codes_df <- ProcWAS::icd10_cpt4_source_to_ccsr   %>% 
select(vocabulary, code)  %>% 
slice_sample(n = 10000)  %>% 
mutate(id = rep(1:1000, each = 10),
       date = seq.Date(Sys.Date(), by = "-1 day", length.out = 10000)
)

# Map source codes to CCSR categories and reshape to wider
phenotypes_df <- create_ccsr_phenotypes(source_codes_df = source_codes_df)

# Create covariates data frame. Some of the individuals in the covariates data frame are not present in phenotypes data frame
covariates_df <- tibble(
  id = 1:15000,
  age = sample(18:65, 15000, replace = TRUE),
  sex = sample(c("M", "F"), 15000, replace = TRUE),
  has_my_variant_of_interest = sample(c(TRUE, FALSE), 15000, replace = TRUE)
)

pheno_cov_df <- merge_pheno_with_cov(phenotypes_df = phenotypes_df,
                                      covariates_df = covariates_df)  %>% 
# create signal for CAR007
mutate(has_my_variant_of_interest = ifelse(CAR007 == TRUE & runif(n()) < 0.4, TRUE, has_my_variant_of_interest))

phewas_df <- apply_sex_specific_exclusions(phewas_df = pheno_cov_df,
                                                   name_of_sex_column = "sex")

results <- phewas_ext(data = phewas_df, 
                     phenotypes = stringr::str_subset(names(phewas_df), "^[A-Z]{3}[0-9]{3}$"), # extract column names with 3 uppercase letters followed by 3 digits consistent with CCSR naming 
                     genotypes=c("has_my_variant_of_interest"),
                     covariates=c("sex", "age"),
                     additive.genotypes = FALSE,
                     cores = 1)

manhattan <- plotManhattan(results,
                           suggestive.line = NA,
                           significant.line = 0.05 / nrow(results),
                           annotate.level = 0.05 / nrow(results),
                           point.size = 2,
                           annotate.size = 3.5)

manhattan
