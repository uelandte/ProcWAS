#' Performs ProcWAS
#' 
#' This function perofrms ProcWAS using generalized linear model with optional adjustment for user-specified covariates
#' 
#' @param phenotypes The names of the outcome variables (i.e., the procedures)
#' @param genotypes The names of the prediction variables (i.e., the genotype)
#' @param data A data frame containing all of the variables to be used in the ProcWAS analysis
#' @param covariates The names of the covariates
#' @param cores The number of cores to be used in the parallel socket cluster implementation.
#' @param additive.genotypes Identifies if additive genotypes are being supplied (default is TRUE).
#' @param return.models Returns a list of the complete models (default is FALSE).
#' @param min.records The minimum number of records necessary to perform a test. For logistic regression, there must be at least this number of cases and at least this number of controls. The default value is 20.
#' @param MASS.confint.level Calculates a confidence interval at the specified level. Logistic models will report OR CIs. Confidence intervals are not calculated by default.
#' @return Returns a data frame containing the results of the PheWAS analysis (e.g., a data frame with the following columns: phenotype, snp, covariates, beta, SE, OR, p, type, n_total, n_cases, n_controls, HWE_p, allele_freq, n_no_snp, formula, expanded_formula, note).
#' @source https://github.com/monikagrabowska/PedsPheWAS/blob/main/R/phewas_ext.R
#' @source https://github.com/PheWAS/PheWAS/blob/master/R/phewas_ext.R
#' @export

phewas_ext <-
  function(phenotypes, genotypes, data, covariates=NA, cores=1, additive.genotypes=T,
           return.models=F, min.records=20, MASS.confint.level=NA) {
    #Require an input data frame
    if(missing(data)) {
      stop("A data frame must be supplied in 'data'")
    }
    #Rename outcomes and predictors parameters if used.
    if(missing(phenotypes)) {
      if(!missing(outcomes)) phenotypes=outcomes
      else stop("Either phenotypes or outcomes must be passed in.")
    }
    if(missing(genotypes)) {
      if(!missing(predictors)) genotypes=predictors
      else stop("Either genotypes or predictors must be passed in.")
    }
    #Convert covariates to a list if it is not one
    if(class(covariates)!="list") { covariates=list(covariates)}
    
    # only method supported for now from original PheWAS documentation is glm
    association_method=phe_as_ext
    
    para=(cores>1)
    
    #Create the list of combinations to iterate over
    full_list=data.frame(t(expand.grid(phenotypes,genotypes,covariates,stringsAsFactors=F)),stringsAsFactors=F)
    
    #If parallel, run the parallel version.
    if(para) {
      #Check to make sure there is no existing phewas cluster.
      if(exists("phewas.cluster.handle")) {
        #If there is, kill it and remove it
        message("Old cluster detected (phewas.cluster.handle), removing...")
        try(stopCluster(phewas.cluster.handle), silent=T)
        rm(phewas.cluster.handle, envir=.GlobalEnv)
      }
      message("Starting cluster...")
      assign("phewas.cluster.handle", makeCluster(cores), envir = .GlobalEnv)
      message("Cluster created, finding associations...")
      clusterExport(phewas.cluster.handle,c("data"), envir=environment())
      clusterCall(phewas.cluster.handle,library,package="dplyr",character.only=T)
      # Loop across every phenotype- iterate in parallel
      result <-parLapplyLB(phewas.cluster.handle, full_list, association_method, additive.genotypes=additive.genotypes,
                           confint.level=MASS.confint.level, min.records=min.records,
                           return.models=return.models)
      # Once we have succeeded, stop the cluster and remove it.
      stopCluster(phewas.cluster.handle)
      rm(phewas.cluster.handle, envir=.GlobalEnv)
    } else {
      # Otherwise, just use lapply.
      message("Finding associations...")
      result=lapply(full_list,FUN=association_method, additive.genotypes=additive.genotypes,
                    confint.level=MASS.confint.level, my.data=data, min.records=min.records,
                    return.models=return.models)
    }
    
    if(return.models) {
      message("Collecting models...")
      models=lapply(result,function(x){attributes(x)$model})
      names(models)=sapply(models,function(x){paste0(as.character(terms(x))[c(2,1,3)],collapse=" ")})
    }
    
    message("Compiling results...")
    successful.phenotypes=na.omit(sapply(result,function(x){attributes(x)$successful.phenotype}))
    n.tests=length(successful.phenotypes)
    successful.phenotypes=unique(successful.phenotypes)
    successful.genotypes=unique(na.omit(sapply(result,function(x){attributes(x)$successful.genotype})))
    sig=bind_rows(result)
    
    #Report warning if any convergence errors
    if(max(grepl(pattern = "[Error: The model did not converge]", sig$note, fixed=TRUE))){
      warning("Not all models converged, check the notes column for details.")
    }
    
    message("Cleaning up...")
    
    if(return.models){sig=list(results=sig,models=models)}
    
    return(sig)
  }