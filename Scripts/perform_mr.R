
######################################################################################################
# This version was updated 17 May 2023
######################################################################################################
#
# "perform_mr" performs two-sample MR analyses using various univariable mendelian randomization methods. It will
#  perform preprocessing including harmonisation, and then applies a selection of MR methods to the 
# data, assesses 
#  heterogeneity across SNPs, and plots a forest plot of causal estimates from the methods selected.
#
#  This function depends on the following R packages:
#  Information on installing the packages is included here but it is easier to use the 'Examples' to set this up!
#
#      + Available in  CRAN: rio, tidyverse,data.table, MendelianRandomization, testthat, gsmr, grid, 
#        gridExtra,ieugwasr,patchwork
#      + Others:
#          + TwoSample MR:
#                 remotes::install_github("MRCIEU/TwoSampleMR")
#          +  gsmr:
#                install.packages("survey")
#                install.packages("markdown")
#
#                install.packages("https://yanglab.westlake.edu.cn/software/gsmr/static/gsmr_1.1.1.tar.gz",repos=NULL,type="source")
#      + Use the function 'p_load' from the R-package pacman to load them all (See example file)
#
#   Available Mendelian Randomization methods:
#         + Assuming independence across SNPs:
#                Inverse variance weighted (fixed effects)
#                Inverse variance weighted (multiplicative random effects)
#                Weighted median
#                Weighted mode
#                MR Egger
#         + Accounting for correlated SNPs
#                Generalised inverse variance weighted (fixed effects)
#                Generalised inverse variance weighted (multiplicative random effects)
#         + Adjusting/accounting for outliers:
#                MR PRESSO adjusted
#                GSMR adjusted
#
#   Arguments of the function:
#       + ss_exposure:      Summary stats for the exposure. 
#                           Required columns: 'rsid', 'effect_allele', 'other_allele', 'beta', 'se', 'pval'. 
#                           However, 'eaf' is needed to subset on MAF and to try to harmonise palindromic SNPs (see harmo_action)
#       + ss_outcome:       Summary stats for the outcome. 
#                           Required columns: 'rsid', 'effect_allele', 'other_allele', 'beta', 'se', 'pval'. 
#                           However, 'eaf' is needed to try to harmonise palindromic SNPs (see harmo_action)
#       + ss_harmo:         A harmonised summary stats file (optional). This will skip preprocessing 
#                           and proceed to analysis. Required columns: 'rsid' (or 'SNP'), 'beta.exposure', 
#                           'se.exposure', 'beta.outcome', 'se.outcome'
#       + cis:              If TRUE, use only cis variants defined by gene_chr, gene_start,gene_end, 
#                           and flank_kb
#       + gene_chr:         Chromosome where the gene of interest is located
#       + gene_start:       Starting position of the gene
#       + gene_end:         Ending position of the gene
#       + flank_kb:         Width of the gene flank (in kilobases)
#       + maf_thresh:       Minor allele frequency (maf) to filter on (default = 0.01)
#       + p_thresh:         P-value threshold to filter on (default = 5x10^-8)
#       + harmo_action      Action for dealing with palindromic SNPs: 1= assume all are on positive 
#                                                                       strand;
#                                                                     2 = try to infer alleles; 
#                                                                     3 = drop palindromic SNPs
#                           See github page for more information.
#       + map_rsids         Map RSID from chr and pos after harmonisation for downstream clumping and LD matrices (see note below). 
#                           Experimental. Default=F
#       + map_build         Build to map RSID (Default = "ch37" or choose "ch38")
#       + clump             If TRUE, clump the variants  (default = TRUE)
#       + clump_kb          Distance for clumping (in kilobases) (default = 250)
#       + clump_r2          r-squared for clumping (default = 0.01) # Note: Based on conversations with FM
#       + clump_pop         population for clumping (default = "EUR"). Options are
#                           'EUR', 'SAS', 'EAS', 'AFR', 'AMR'
#       + clump_only:       If TRUE, clump the variants BEFORE harmonising and just produce the list (default= FALSE)
#       + name_exp          Name of the exposure (optional) (default = 'Exposure')
#       + name_outcome      Name of the outcome (optional) (default = 'Outcome')
#       + which_methods     Which MR methods to perform
#                           Default =  'all' to run all methods;
#                                      'basic' to run ivw_fe, ivw_mre, w_median, w_mode, mr_egger;
#                                       or specify among:
#
#                                mr_wald_ratio: Wald ratio (single variant analyses)
#                                mr_ivw_fe: IVW with fixed effects
#                                mr_ivw_mre: IVW with multiplicative random effects
#                                mr_ivw_corr_fe: Generalised IVW for correlated variants with fixed effects
#                                mr_ivw_corr_mre: Generalised IVW for correlated variants with multiplicative 
#                                              random effects
#                                mr_weighted_median: weighted median
#                                mr_weighted_mode: weighted mode
#                                mr_egger_regression:  MR-Egger 
#                                mr_presso_unadj: MR PRESSO without outlier correction
#                                mr_presso_adj: MR PRESSO with outlier correction
#                                gsmr_unadj: GSMR without removing Heidi outliers
#                                gsmr_adj: GSMR after removing Heidi outliers
#                                The function will always include ivw_fe or ivw_mre regardless of selection!
#          
#       + ld_mat:            LD matrix required for generalised IVW or GSMR.
#                            Default = user supplied matrix with no missing values.
#                            The matrix should have rows/columns in the same 
#                            order as the exposure summary statistics and the column names should
#                            correspond to the rsids. See https://ibb.co/6Bb5hJc for an example.
#                            Or use "try_gen' to try to generate a matrix (for up to 500 
#                            variants)
#       + ld_mat_pop         Population from 1000G used to generate the matrix (default = 'EUR').Options are
#                           'EUR', 'SAS', 'EAS', 'AFR', 'AMR'
#       + presso_nboot:      Number of bootstraps for MR PRESSO
#                            Default = 'min_rec' to calculate and use the minimum recommended 
#                             by the package. Or specify a real number (e.g. 500)
#       + which_plots:       MR methods to include in forest plot
#                            Default = "all" to plot all methods;
#                                      'none' to include no plot;
#                                       or specify among methods in selected in 'which_methods'
#
#       + outcome_binary: Flag to identify whether the outcome is binary (T) or continuous (F)
#                         (for presentation of results)
#
#   Outputs:
#       + An R list with 6 objects
#           + results:  Data frame with causal estimates by method
#           + het:      Data frame with results from heterogeneity tests (Mean F-statistic, Q-statistic, 
#                       Rucker's Q difference and MR-PRESSO global test)
#           + exclude:  Data frame with rsids of SNPs excluded by MR-PRESSO and GSMR (if relevant)
#           + for_plot: ggplot object with forest plot
#           + harmo_df: Data frame with the harmonised data
#           + harmo_exclude: Data frame with list of variants excluded during harmonisation and
#             the reason
#
# Note: The forest plot is provided as a rough visual aid and is not designed
#       to accommodate every use case. Please use the results table provided to create your own
#       plots.
#
# Please read the ReadMe file at https://github.com/adamkvonende/MRpipeline!
#
# 
######################################################################################################


# Function
perform_mr <-function(ss_exposure, ss_outcome,ss_harmo, cis = F, gene_chr, gene_start,gene_end,flank_kb=500,
                         maf_thresh=0.01, p_thresh=5*10^-8, 
                          harmo_action=3, map_rsids=F, map_build="ch37",
                         clump=TRUE, clump_kb=250, clump_r2 = 0.01,clump_pop="EUR", clump_only=F,
                          name_exp="Exposure", name_outcome="Outcome",
                          which_methods="all", ld_mat ="try_gen", ld_mat_pop="EUR",
                          presso_nboot = "min_rec",which_plots="all", 
                          outcome_binary=T) {
       
  # Switch for whether harmonised file is include or not
  if (!missing(ss_harmo)) {
    data_harm<-ss_harmo
    req.columns.harmo <- c("rsid","beta.exposure","se.exposure","beta.outcome","se.outcome")

    if (!all(req.columns.harmo %in% names(data_harm)) && !("SNP" %in%  names(data_harm) )) {
      stop("The following columns in the harmonised file are required: rsid, beta.exposure, se.exposure, beta.outcome, se.outcome. Stopping...\n\n")
    }
    
    
    if (!"SNP" %in% names(data_harm)& "rsid" %in% names(data_harm)) {data_harm$SNP=data_harm$rsid}
    if (!"id.exposure" %in% names(data_harm)) {data_harm$id.exposure=name_exp}
    if (!"id.outcome" %in% names(data_harm)) {data_harm$id.outcome=name_outcome}
    if (!"exposure" %in% names(data_harm)) {data_harm$exposure=name_exp}
    if (!"outcome" %in% names(data_harm)) {data_harm$outcome=name_outcome}
    if (!"mr_keep" %in% names(data_harm)) {data_harm$mr_keep=T}
    cat("A harmonised file is provided so preprocessing will not be performed","\n\n")
    
  } else {
    cat("*** The function is now preprocessing the data for analysis ***\n\n")
    
    req.columns <- c("rsid","effect_allele","other_allele","beta","se", "pval")
    non.req.columns<- c("chr","pos", "eaf")
    
    if (any(!req.columns %in% names(ss_exposure))) {
      stop("The following columns in the exposure summary stats are required: rsid, effect_allele, other_allele, beta, se, pval. Stopping...")
    }
    data_exp <- ss_exposure 
    for (col_name in non.req.columns) {
      if (!(col_name %in% names(data_exp))) {
        data_exp[[col_name]] <- NA
      }
    }
    
    data_exp <- data_exp %>% dplyr::mutate_at(c("chr","effect_allele","other_allele"), function(x) as.character(x))
    data_exp <- data_exp %>% dplyr::mutate_at(c("pos", "beta", "se", "pval", "eaf"), function(x) as.numeric(x))
    data_exp2<- data_exp %>% dplyr::filter(pval < p_thresh) # Subset on p
    cat(paste(nrow(data_exp2), "variants were retained in the exposure data after removing variants with p-value >=",p_thresh,"\n\n" ))
    if (sum(is.na(data_exp2$eaf)) >0) {
    cat(paste(sum(is.na(data_exp2$eaf)), "variants have missing values for EAF and will not be filtered on maf \n\n"))
    data_exp2$eaf[is.na(data_exp2$eaf)]<-.5
    }
    data_exp3 <- data_exp2 %>% dplyr::filter(eaf>=maf_thresh&eaf<=1-maf_thresh) # Subset to maf 
    cat(paste(nrow(data_exp3), "variants were retained in the exposure data after removing variants with maf <", maf_thresh,"\n\n"))
    data_exp4 <- data_exp3  %>% filter(!(effect_allele==other_allele)) %>% filter(!(nchar(effect_allele)==0|nchar(other_allele)==0))
    cat(paste(nrow(data_exp4), "variants were retained in the exposure data after removing variants with missing or identical alleles","\n\n"))
    
    
      if (cis ==T) {
     data_exp4 <- data_exp4  %>%  mutate(cis =ifelse(as.character(gene_chr)==chr & pos >= gene_start-flank_kb*1000 & pos<=gene_end+flank_kb*1000,1,0 ),
                                         cis = ifelse(is.na(cis),0, cis))
     if ((sum(data_exp4$cis)==0)) {
       stop("No cis variants found: stopping...") 
       } else {
        cat(paste(sum(data_exp4$cis), "cis variants were found in the exposure data and carried forward for analysis","\n\n" ))
         data_exp4 <- data_exp4 %>% filter(cis==1)
       } 
      } 
    
    # Convert data sets to TwoSampleMR versions
    exp_in <- data_exp4  %>% transmute(SNP=rsid, beta.exposure=beta, se.exposure=se, pval.exposure=pval,effect_allele.exposure=effect_allele, other_allele.exposure=other_allele, eaf.exposure=eaf,exposure=name_exp,id.exposure="")
    
    
    # Clumping only
    
    
    
    if (clump_only==T) {
      cat("*** The function is now performing clumping on the pre-harmonised data; this will stop when complete!***\n\n")
      
      # Check for server side problems
      clump.capture <- capture_messages(clump_data(exp_in,clump_kb = clump_kb,clump_r2 = clump_r2,clump_p1 = 1, clump_p2 = 1,pop = clump_pop))
      fail.capture <- grep("Failed", clump.capture, value = TRUE)
      if (length(fail.capture)>0) {
        stop("MRbase failed to retrieve results from the server. Clumping failed and will return 0 clumped variants. Please try again.")
      }
      # Perform clumping
      clump_out <- suppressMessages(clump_data(exp_in,
                                            clump_kb = clump_kb,
                                            clump_r2 = clump_r2,
                                            clump_p1 = 1, clump_p2 = 1,
                                            pop = clump_pop))
      clumped_list <- clump_out %>% dplyr::select(rsid=SNP)
      
      cat("The function will now stop and return ONLY a list of clumped variants.\n")
      cat(paste(nrow(clump_out ), "variants were output after clumping at r2 =", clump_r2,"\n\n"))
      
      return(clumped_list)

      
    } 
    
    
     
    if (any(!req.columns %in% names(ss_outcome))) {
      stop("The following columns in the outcome summary stats are required: rsid, effect_allele, other_allele, beta, se, pval. Stopping...")
    }
    data_out <- ss_outcome 
    for (col_name in non.req.columns) {
      if (!(col_name %in% names(data_out))) {
        data_out[[col_name]] <- NA
      }
    }
    data_out <- data_out %>% dplyr::mutate_at(c("chr","effect_allele","other_allele"), function(x) as.character(x))
    data_out <- data_out %>% dplyr::mutate_at(c("pos", "beta", "se", "pval", "eaf"), function(x) as.numeric(x))
    # Convert data sets to TwoSampleMR version
    outcome_in <- data_out %>%  transmute(SNP=rsid, beta.outcome=beta, se.outcome=se, pval.outcome=pval,effect_allele.outcome=effect_allele, other_allele.outcome=other_allele, eaf.outcome=eaf,outcome=name_outcome,id.outcome="")
    
    
  
    ### Harmonise data
  

    cat("*** The function is now performing harmonisation ***\n\n")
    
    
        if (!is.na(sum(exp_in$SNP == ""))&sum(exp_in$SNP=="")>0) {
      cat(paste(sum(exp_in$SNP==""), "variants in the exposure data set have missing rsids and will be dropped before harmomisation","\n\n"))
      exp_in <- exp_in %>% filter(SNP!="")
      }
    
    not_in_outcome_ss <-  exp_in$SNP[!exp_in$SNP %in% outcome_in$SNP]
    if (!is_empty(not_in_outcome_ss)) {
    cat(length(not_in_outcome_ss), "SNPs were not found in the outcome data set according to rsid. Consider finding proxies.See 'harmo_exclude' for the list.\n\n")
    df.exclude.outcome = data.frame(SNP = not_in_outcome_ss, reason = "Not in outcome data set")
      }

    data_harm<- suppressMessages(harmonise_data(exp_in, outcome_in, action = harmo_action)) #'action:1=assume all alleles on forward strand; 2=try to infer alleles;3=drop palindromic
    if (nrow(data_harm)==0) {stop("There are no variants remaining after attempting to harmonise")}
    
    if (harmo_action==1) {
      cat(paste("All variants were assumed to be on the positive strand. There were",  sum(data_harm$palindromic), 
                "palindromic variants and none were dropped","\n\n"))
    }
    if (harmo_action==2) {
      cat(paste("There were",  sum(data_harm$palindromic), "palindromic variants and",
                  sum(data_harm$palindromic)-sum(data_harm$ambiguous), "were retained after excluding ambiguous variants","\n\n" ))
    }
    
    if (harmo_action==3) {
      cat("There were",  sum(data_harm$palindromic), "palindromic variants; these were dropped\n\n")
    }
    
    # Get list of SNPs with incompatible alleles
    harm.capture <- capture_messages(harmonise_data(exp_in, outcome_in, action = harmo_action))
    incompat.capture <- grep("compatible", harm.capture, value = TRUE)
    if (!is_empty(incompat.capture)) {
      incompat.capture2  <- strsplit(incompat.capture, ":")[[1]][[2]]
      incompat.capture3  <- strsplit(incompat.capture2, ", ")
      incompat.snps <-incompat.capture3[[1]]
      incompat.snps <- gsub("\n", "", incompat.snps)
      cat(length(incompat.snps), "SNPs were excluded due to incompatible alleles. See 'harmo_exclude' for the list")
      df.exclude.incompat = data.frame(SNP = incompat.snps, reason = "Incompatible alleles")
    }
    
    harmo_exclude <- bind_rows(get0("df.exclude.outcome"), get0("df.exclude.incompat"))

    # Filter to those to keep
    
    data_harm <- data_harm %>% filter(mr_keep)
    cat("\n\nThere are", nrow(data_harm), "variants available for analysis after harmonisation\n\n")
    if (nrow(data_harm)==0) {
      stop("No variants remaining after harmonisation; stopping")
    }
    
 
    # Map rsids 
    
    if (map_rsids==T) {
      tmp1 <- tidyr::separate(data_harm, SNP, into = c("chr", "pos"), sep = ":")
      tmp1  <- Map_rsid_from_chrpos(tmp1 ,build=map_build)
      tmp1  <- tmp1 %>% dplyr::mutate(SNP = rsid_mapped)
      data_harm <- tmp1  %>% dplyr::select(SNP, all_of(c("effect_allele.exposure", "other_allele.exposure", 
                                                     "effect_allele.outcome", "other_allele.outcome", "beta.exposure", 
                                                     "beta.outcome", "eaf.exposure", "eaf.outcome", "remove", "palindromic", 
                                                     "ambiguous", "id.outcome", "se.outcome", "pval.outcome", "outcome", 
                                                     "se.exposure", "pval.exposure", "exposure", "id.exposure", "action", 
                                                     "mr_keep", "samplesize.outcome")))

      if (nrow(data_harm)==0) {
        stop("No variants remaining; stopping")
      }
    }
    
    
    
    
    ### Clumping after harmo
    
    if (clump==T & clump_only !=TRUE) {
  
          cat("*** The function is now performing clumping on the harmonised data***\n\n")
          # Check for server side problems
          clump.capture <- capture_messages(clump_data(exp_in,clump_kb = clump_kb,clump_r2 = clump_r2,clump_p1 = 1, clump_p2 = 1,pop = clump_pop))
          fail.capture <- grep("Failed", clump.capture, value = TRUE)
          if (length(fail.capture)>0) {
            stop("MRbase failed to retrieve results from the server. Clumping failed and will return 0 clumped variants. Please try again.")
          }
          data_harm <- suppressMessages(clump_data(data_harm,
                                  clump_kb = clump_kb,
                                  clump_r2 = clump_r2,
                                  clump_p1 = 1, clump_p2 = 1,
                                  pop = clump_pop))
          
          cat(paste(nrow(data_harm), "variants were retained after clumping at r2 =", clump_r2,"\n\n"))
          if (nrow(data_harm)==0) {
            stop("No variants available after clumping; stopping")
          }
      
      
                      } 
                                
  
  }

  cat("*** The function is now performing MR on the harmonised data *** \n\n")
  
  
  if (sum(which_methods=="basic")==1) {
    which_methods<-  c("mr_ivw_fe","mr_ivw_mre",
                       "mr_weighted_median","mr_weighted_mode","mr_egger_regression")
  } else if (sum(which_methods=="all")==1) {
    which_methods<-  c("mr_ivw_fe","mr_ivw_mre","mr_ivw_corr_fe","mr_ivw_corr_mre",
                       "mr_weighted_median","mr_weighted_mode","mr_egger_regression","mr_presso_unadj",
                       "mr_presso_adj","gsmr_unadj","gsmr_adj")
  } else {
    which_methods<-unique(c("mr_ivw_fe","mr_ivw_mre",which_methods)) # Always include IVW 
  }
  
  
  # If SNPs = 1 then only run Wald ratio
  
  if (nrow(data_harm)==1) {
    cat("There is only 1 variant; only the Wald ratio method will be performed!","\n\n")
    cat("Generalised IVW requires > 1 variants; skipping this method","\n")
    #which_methods<-setdiff(which_methods, c("mr_ivw_corr_fe","mr_ivw_corr_mre"))
    which_methods<-"mr_wald_ratio"
    
  }
  
  
  # If SNPs <=2 ignore Egger, median, and mode
  
  if (nrow(data_harm)<=2) {
    cat("Weighted median requires > 2 variants; skipping this method\n")
    cat("Weighted mode requires > 2 variants; skipping this method\n")
    cat("MR Egger requires > 2 variants; skipping this method\n")
    which_methods<-setdiff(which_methods, c("mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
    
  }
  
  
  # If SNPs <=3 ignore PRESSO
  
  if (nrow(data_harm)<=3) {
    cat("MR PRESSO requires > 3 variants; skipping this method","\n")
    which_methods<-setdiff(which_methods, c("mr_presso_unadj","mr_presso_adj"))
    
  }
  
  # If SNPs <= 10 ignore GSMR
  
  if (nrow(data_harm)<10) {
    cat("GSMR requires >= 10 variants; skipping this method","\n")
    which_methods<-setdiff(which_methods, c("gsmr_unadj","gsmr_adj"))
    
  } 
  
  # If plot selections not a subset of method selections
  
  
  if (sum(which_plots!="none")>=1&sum(which_plots!="all")>=1& !all(which_plots %in% which_methods )) {
    stop("You have selected results to be plotted that are not included in the results. 
          This could be because the method could not be run given the data you have provided. 
          Please change your plot selections.")
  }
  
  
  

  
  if (any(c("mr_ivw_corr_fe", "mr_ivw_corr_mre","gsmr_unadj","gsmr_adj")  %in% which_methods)) {
    
    if  (sum(ld_mat == "try_gen")>=1) {
      cat("\nThe function is attempting to generate an LD matrix for the harmonised data...","\n\n")
      ld_mat <- suppressMessages(ieugwasr::ld_matrix(variants=data_harm$SNP,pop = ld_mat_pop))
        if (dim(ld_mat)[1] < length(data_harm$SNP)) {
          stop("At least one SNP was not found in the reference panel; stopping. Please submit your own LD matrix if you want to run Generalised IVW or GSMR.")
        }   else {cat("The LD matrix has been successfully created.\n\n") }

    }
            else {
      ld_mat <-ldmat
      cat("Using the LD matrix supplied by the user","\n\n")
      
    }
    
    if (missing(ld_mat)) {
      stop("You must supply an LD matrix; stopping")
    }
    
    if (any(is.na(ld_mat))) {
      stop("The supplied LD matrix has missing values; stopping")
    }
    
    # Define number of individuals
    
    ld_mat_pop_values <- c('AFR', 'AMR','EAS', 'SAS', 'EUR')
    ld_mat_n_values <- c(661, 347, 504, 489, 503)
    ld_mat_n <- ld_mat_n_values[match(ld_mat_pop, ld_mat_pop_values)]
    

  }
    

  
  ### Basic MR methods
  base_methods <- c("mr_wald_ratio","mr_ivw_fe","mr_ivw_mre","mr_weighted_median","mr_weighted_mode", "mr_egger_regression")
  methods_to_use <-  base_methods[base_methods %in% which_methods]
  
  
  if (nrow(data_harm)==1) {
    methods_to_use <-"mr_wald_ratio"
  }
  
  
  suppressMessages(mr.results <- mr(data_harm, method_list=methods_to_use) %>% dplyr::select(method, nsnp, b, se, pval))

  
 
  
  ### Generalised IVW 


  if ("mr_ivw_corr_fe" %in% which_methods) {
    
    mrinput <- with(data_harm,mr_input(bx=beta.exposure,bxse=se.exposure,
                                  by=beta.outcome,byse=se.outcome,
                                  correlation = as.matrix(ld_mat),
                                  outcome=outcome,exposure = exposure,
                                  snps=SNP,
                                  effect_allele = effect_allele.exposure,
                                  other_allele = other_allele.exposure))
    
    
    corrivw_fe.out <-  tryCatch({
      MendelianRandomization::mr_ivw(mrinput,method="fixed")},
      error = function(e) paste("Error: ", e))  
    
  
    if (is.character(corrivw_fe.out)) {
      cat("mr_ivw_corr_fe produced an error; skipping\n\n")
    } else {
    
    corrivw_fe.results <- data.frame(
      method="Generalised inverse variance weighted (fixed effects)",
      nsnp=corrivw_fe.out@SNPs,
      b=corrivw_fe.out@Estimate,
      se=corrivw_fe.out@StdError,
      pval=corrivw_fe.out@Pvalue)
    }
    
  }
  
  if ("mr_ivw_corr_mre" %in% which_methods) {
    
    mrinput <- with(data_harm,mr_input(bx=beta.exposure,bxse=se.exposure,
                                  by=beta.outcome,byse=se.outcome,
                                  correlation = as.matrix(ld_mat),
                                  outcome=outcome,exposure = exposure,
                                  snps=SNP,
                                  effect_allele = effect_allele.exposure,
                                  other_allele = other_allele.exposure))

    
    corrivw_mre.out <-  tryCatch({
      MendelianRandomization::mr_ivw(mrinput,method="fixed")},
      error = function(e) paste("Error: ", e))  
    
    
    if (is.character(corrivw_mre.out)) {
      cat("mr_ivw_corr_mre produced an error; skipping\n\n")
    } else {

    corrivw_mre.results <- data.frame(method="Generalised inverse variance weighted (multiplicative random effects)",
                                      nsnp=corrivw_mre.out@SNPs,
                                      b=corrivw_mre.out@Estimate,
                                      se=corrivw_mre.out@StdError,
                                      pval=corrivw_mre.out@Pvalue)
    }
  }
  
  
  ### MR PRESSO

  
  if ("mr_presso_unadj" %in% which_methods |"mr_presso_adj" %in% which_methods) {
    
    
    if (presso_nboot == "min_rec") {
      min_nb <- ceiling(nrow(data_harm)/0.05) 
      cat(paste("Running MR-PRESSO using", min_nb, "bootstraps, the minimum recommended for a stable outlier test according to documentation" ,"\n\n" ))
      
    } else {
      min_nb = presso_nboot
      cat(paste("Running MR-PRESSO using", min_nb, "bootstraps specified by the user.","\n\n"  ))
          if (presso_nboot < ceiling(nrow(data_harm)/0.05)) {
            message("The number of bootstraps requested for the PRESSO outlier test is less than recommended and may produce an error. ")
          }
      
    }
    

    
    mrpresso.out <- run_mr_presso(data_harm,NbDistribution = min_nb)
    mrpresso_global_p <- as.character(mrpresso.out [[1]][["MR-PRESSO results"]][["Global Test"]][["Pvalue"]])
    
  }
  

  if ("mr_presso_unadj" %in% which_methods) {
    
    mrpresso_unadj.results <- mrpresso.out[[1]][["Main MR results"]][1,]
    colnames(mrpresso_unadj.results) <- c("exposure","MR analysis","b","se","t-stat","pval")
    mrpresso_unadj.results<- mrpresso_unadj.results %>% transmute(method = "MR PRESSO unadjusted",nsnp=nrow(data_harm),b,se, pval )
    
    
  }
  
  if ("mr_presso_adj" %in% which_methods) {
    
    mrpresso_adj.results <- mrpresso.out[[1]][["Main MR results"]][2,]
    
    #if no outliers then estimates will be NA
    #take results from raw model instead
    if (is.na(mrpresso_adj.results$`Causal Estimate`)) {
      mrpresso_adj.results <-  mrpresso.out[[1]][["Main MR results"]][1,]
    }
    
    colnames(mrpresso_adj.results) <- c("exposure","MR analysis","b","se","t-stat","pval")
    #Calculate excluded SNPs
    excluded_presso_index <- mrpresso.out[[1]][["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]
    if (!is.null(excluded_presso_index)) {
      exclude1 <- data.frame(SNP=data_harm$SNP[excluded_presso_index],excluded="PRESSO")
    }
    mrpresso_adj.results<- mrpresso_adj.results %>% transmute(method = "MR PRESSO adjusted",nsnp=length(data_harm$SNP) - length(excluded_presso_index),b,se, pval )
    
    
    
  }
  
  ### GSMR

  
  if ("gsmr_unadj" %in% which_methods) {
    
    rownames(ld_mat) <- colnames(ld_mat)
    suppressMessages(gsmr_unadj.out <- with(data_harm, gsmr(bzx=beta.exposure,bzx_se=se.exposure,bzx_pval=pval.exposure,
                                                       bzy=beta.outcome,bzy_se=se.outcome,bzy_pval = pval.outcome,
                                                       ldrho=as.matrix(ld_mat),snpid=colnames(ld_mat),n_ref = ld_mat_n,
                                                       heidi_outlier_flag = F,gwas_thresh = 1,ld_r2_thresh=1,nsnps_thresh=1)))
    
    
    gsmr_unadj.results <- data.frame(method="GSMR unadjusted",
                                     nsnp=length(gsmr_unadj.out$used_index),
                                     b=gsmr_unadj.out$bxy,
                                     se=gsmr_unadj.out$bxy_se,
                                     pval=gsmr_unadj.out$bxy_pval)
    
  }
  
  if ("gsmr_adj" %in% which_methods) {
    
    rownames(ld_mat) <- colnames(ld_mat)
    suppressMessages(gsmr_adj.out <- with(data_harm, gsmr(bzx=beta.exposure,bzx_se=se.exposure,bzx_pval=pval.exposure,
                                                     bzy=beta.outcome,bzy_se=se.outcome,bzy_pval = pval.outcome,
                                                     ldrho=as.matrix(ld_mat),snpid=colnames(ld_mat),n_ref = ld_mat_n,
                                                     gwas_thresh = 1,ld_r2_thresh=1,nsnps_thresh=1)))
    
    excluded_snps_gsmr <- data_harm$SNP[-gsmr_adj.out[["used_index"]]]
    if (!rlang::is_empty(excluded_snps_gsmr)) {
      exclude2 <- data.frame(SNP=excluded_snps_gsmr,excluded="GSMR")
    }
    gsmr_adj.results <- data.frame(method="GSMR adjusted",
                                   nsnp=length(gsmr_adj.out$used_index),
                                   b=gsmr_adj.out$bxy,
                                   se=gsmr_adj.out$bxy_se,
                                   pval=gsmr_adj.out$bxy_pval)
    
  }
  

  #Combine different results tables
  df.out <- bind_rows(get0("mr.results"), get0("corrivw_mre.results"),get0("corrivw_fe.results"),
                      get0("mrpresso_unadj.results"), get0("mrpresso_adj.results"),
                      get0("gsmr_unadj.results"), get0("gsmr_adj.results"))
  
  
  if (outcome_binary==T) {
    
    results.table  <- df.out  %>% transmute(outcome=data_harm$outcome[1], 
                                            exposure= data_harm$exposure[1],
                                            method,nsnp,b,se,
                                            or = exp(b), 
                                            ll=exp(b-qnorm(0.975)*se),
                                            ul=exp(b+qnorm(0.975)*se),
                                            pval)
    
  } else {
    results.table  <- df.out  %>% transmute(outcome=data_harm$outcome[1], 
                                            exposure= data_harm$exposure[1],
                                            method,nsnp,b,se,
                                            ll=b-qnorm(0.975)*se,
                                            ul=b+qnorm(0.975)*se,
                                            pval)
    
    
  }

  
  ### Create het stats output table

  
  beta_format <- function(x){if (x<0.01) {formatC(x,format="e",digits=2)}
    else{ifelse(is.na(x),"",formatC(x,format="f",digits=2))}}
  p_format <- function(x){if (is.na(x)==F & x<0.001) {formatC(x,format="e",digits=2)}
    else{ifelse(is.na(x),"",formatC(x,format="f",digits=3))}}
  
  
  # Mean F-statistic
  mean_fstat = mean((data_harm$beta.exposure^2)/(data_harm$se.exposure^2))
  het1 <- data.frame(Statistic="Mean F-statistic", DF=NA, Estimate=beta_format(mean_fstat) , Pvalue="")
  
  
  # Q-statistic
  
  if (nrow(data_harm)<2) {
    cat("The Q-statistic based heterogeneity test requires > 1 variants; this will not be present in the het table","\n")
  } else {
    het_test <- mr_heterogeneity(data_harm, method_list = c("mr_ivw"))
    het2 <- data.frame(Statistic="Q statistic", DF=het_test$Q_df, Estimate= beta_format(het_test$Q), Pvalue=p_format(het_test$Q_pval))
  }
  
  
  # Rucker's Q difference
  
  if (nrow(data_harm)<3) {
    cat("Rucker's Q requires >= 3 variants; this will not be present in the het table","\n")
  } else {
    rucker_qstats <- mr_rucker(data_harm)[[1]][["Q"]]
    het3 <- data.frame(Statistic="Rucker's Q", DF=rucker_qstats[3,3], Estimate= beta_format(rucker_qstats[3,2]), Pvalue=p_format(rucker_qstats[3,4]))
  }
  
  # Egger intercept
  
  
  if (nrow(data_harm)<3) {
    cat("The Egger intercept requires >= 3 variants; this will not be present in the het table","\n")
  } else {
    x <- mr_pleiotropy_test(data_harm)
    egger_int <- mr_pleiotropy_test(data_harm)$egger_intercept
    egger_p <- mr_pleiotropy_test(data_harm)$pval
    het4 <- data.frame(Statistic="Egger intercept", DF=NA, Estimate= beta_format(egger_int), Pvalue=  p_format(egger_p))
  }
  
  # MR-PRESSO global test
  
  if (exists("mrpresso_global_p")) {
    het5 <- data.frame(Statistic="PRESSO global test", DF=NA, Estimate= NA, Pvalue=  mrpresso_global_p)
  }
  
  # Combine #
  het.table <- bind_rows(get0("het1"), get0("het2"),get0("het3"),get0("het4"), get0("het5"))
  
  
  
  ### Create excluded SNPs output table

  
  # if (exists("excluded_presso_index")&!is.null(excluded_presso_index)) {
  #   exclude1 <- data.frame(SNP=data$SNP[excluded_presso_index],excluded="PRESSO")
  # }
  
  
  exclude.table <- bind_rows(get0("exclude1"), get0("exclude2")) 
  
  
  
  ### Construct forest plot


  
  if(sum(which_plots!="none")>=1){
    
    table <- results.table
    match_methods.df <- data.frame(method_short=c("mr_wald_ratio","mr_ivw_fe","mr_ivw_mre","mr_ivw_corr_fe","mr_ivw_corr_mre","mr_weighted_median","mr_weighted_mode","mr_egger_regression","mr_presso_unadj","mr_presso_adj","gsmr_unadj","gsmr_adj"),
                                   method=c("Wald ratio", "Inverse variance weighted (fixed effects)", "Inverse variance weighted (multiplicative random effects)", "Generalised inverse variance weighted (fixed effects)","Generalised inverse variance weighted (multiplicative random effects)","Weighted median", "Weighted mode", "MR Egger","MR PRESSO unadjusted", "MR PRESSO adjusted", "GSMR unadjusted", "GSMR adjusted"),
                                   method_label = c("Wald ratio", "IVW (fixed effects)", "IVW (random effects)", "Generalised IVW (fixed effects)","Generalised IVW (random effects)","Weighted median", "Weighted mode", "MR Egger","MR PRESSO unadj", "MR PRESSO adj", "GSMR unadj", "GSMR adj")
    )
    
    
    table <- table %>% left_join(match_methods.df, by="method")
    
    
    
    if(outcome_binary==T){table$lab <- paste0(formatC(table$or,digits=2, format="f")," (",formatC(table$ll,digits=2, format="f"),", ",formatC(table$ul,digits=2, format="f"),")")
    } else {table$lab <- paste0(formatC(table$b,digits=2, format="f")," (",formatC(table$b-qnorm(0.975)*table$se,digits=2, format="f"),", ",formatC(table$b+qnorm(0.975)*table$se,digits=2, format="f"),")")}
    
    
    
    table$Index <- 1:nrow(table)
    table$wt <- 1/table$se
    table$plab <- sapply(table$pval,p_format)
    exposure <- toupper(table$exposure[1])
    outcome <- toupper(table$outcome[1])
    
    #which methods to plot
    if(sum(which_plots!="all")>=1){
      table <- subset(table, method_short %in% which_plots)
    } else {
      table <- subset(table, method_short %in% which_methods)
    }
    
    
    
    if(outcome_binary==T) {
      
      #get limits for x-axis
      xmin <- min(table$ll)
      xmin <- round(xmin-0.1,1)
      
      xmax <- max(table$ul)
      xmax <- round(xmax+0.1,1)
      
      p <- ggplot(data=table, aes(y=Index, x=rev(or), xmin=rev(ll), xmax=rev(ul)))+ 
        
        #this adds the ORs to the plot
        geom_point(shape=15, aes(size=rev(wt), y=Index)) +
        #scale_size_continuous(range = c(2,3)) +
        
        #adds the CIs
        geom_errorbarh(height=.1,aes(y=Index)) +
        
        #adding a vertical line at the OR=1 mark
        geom_vline(xintercept=1, color="black", alpha=.5)+
        
        scale_x_continuous(limits=c(xmin,xmax),trans="log10") + 
        
        theme_classic()+
        theme(text=element_text(size=12, color="black"))+
        theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(color="black")) +
        theme(legend.position = 'none',plot.margin=unit(c(5,0,0,0),"pt")) +
        xlab(paste0("Odds ratio per unit higher")) +
        ylab("") +
        ggtitle("\n")
      
      ###################################
      #adding columns of data
      table_base <- ggplot(data=table, aes(y=Index, x=rev(or)))+ 
        ylab(NULL) + xlab("  ") + 
        theme(plot.title = element_text(hjust = 0.5,face="bold"), 
              axis.text.x = element_text(color="white",hjust=0, size = 32), ## This is used to help with alignment
              axis.line = element_blank(),
              axis.text.y = element_blank(), 
              axis.ticks = element_blank(),
              axis.title.y = element_blank(), 
              legend.position = "none",
              panel.background = element_blank(), 
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              
              panel.grid.minor = element_blank(), 
              plot.background = element_blank(),
              plot.margin=unit(c(10,0,0,0),"pt"))
      
      #adding acceleration quintiles and model names
      tab1 <- table_base + 
        ggtitle("Method\n") +
        geom_text(aes(y = Index, x = 0, label = rev(method_label)),hjust=0) + #left align column
        xlim(0,1) +
        theme(plot.title = element_text(hjust=0.07)
        ) #title directly over labels
      
      tab2 <- table_base + 
        ggtitle("OR (95% CI)\n") +
        geom_text(aes(y = rev(Index), x = 0, label = lab))
      
      tab3 <- table_base + 
        ggtitle("P-value\n") +
        geom_text(aes(y = rev(Index), x = 0, label = plab))
      
      
    } else { #if outcome continous
      
      #x-axis limits
      xmin <- min(table$b-qnorm(0.975)*table$se)
      xmin <- round(xmin-0.1,1)
      
      xmax <- max(table$b+qnorm(0.975)*table$se)
      xmax <- round(xmax+0.1,1)
      
      p <- ggplot(data=table, aes(y=Index, x=rev(b), xmin=rev(b-qnorm(0.975)*se), xmax=rev(b+qnorm(0.975)*se)))+ 
        
        #this adds the ORs to the plot
        geom_point(shape=15, aes(size=rev(wt), y=Index)) +
        #scale_size_continuous(range = c(2,3)) +
        
        #adds the CIs
        geom_errorbarh(height=.1,aes(y=Index)) +
        
        #adding a vertical line at the beta=0 mark
        geom_vline(xintercept=0, color="black", alpha=.5)+
        
        theme_classic()+
        theme(text=element_text(size=12, color="black"))+
        theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(color="black")) +
        theme(legend.position = 'none',plot.margin=unit(c(5,0,0,0),"pt")) +
        xlab(paste0("Beta coefficient per unit higher")) +
        ylab("") +
        ggtitle("\n")
      
      
      ###################################
      #adding columns of data
      table_base <- ggplot(data=table, aes(y=Index, x=rev(b)))+ 
        ylab(NULL) + xlab("  ") + 
        theme(plot.title = element_text(hjust = 0.5,face="bold"), 
              axis.text.x = element_text(color="white",hjust=0, size = 32), ## This is used to help with alignment
              axis.line = element_blank(),
              axis.text.y = element_blank(), 
              axis.ticks = element_blank(),
              axis.title.y = element_blank(), 
              legend.position = "none",
              panel.background = element_blank(), 
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              
              panel.grid.minor = element_blank(), 
              plot.background = element_blank(),
              plot.margin=unit(c(10,0,0,0),"pt"))
      
      tab1 <- table_base + 
        ggtitle("Method\n") +
        geom_text(aes(y = Index, x = 0, label = rev(method_label)),hjust=0) + #left align column
        xlim(0,1) +
        theme(plot.title = element_text(hjust=0.07)) #title directly over labels
      
      tab2 <- table_base + 
        ggtitle("Beta (95% CI)\n") +
        geom_text(aes(y = rev(Index), x = 0, label = lab))
      
      tab3 <- table_base + 
        ggtitle("P-value\n") +
        geom_text(aes(y = rev(Index), x = 0, label = plab))
      
    }
    
    #library(patchwork)
    tab1+p+tab2+tab3+plot_layout(nrow=1)
    
    #put everything together
    for_plot <- grid.arrange(tab1,p,tab2,tab3,
                             nrow=1,
                             top=textGrob(paste0("Effect of ",outcome," per unit higher ",exposure),gp=gpar(fontsize=15)))
    
    
  }
  cat("*** The function has finished performing MR for", name_exp, "and", name_outcome, "***")

  
  if(sum(which_plots!="none")>=1) {
    return(list(results=results.table,het=het.table,exclude=exclude.table, harmo_df=data_harm,harmo_exclude=harmo_exclude, ld_mat = ld_mat,for_plot=for_plot))
  } else {
    return(list(results=results.table,het=het.table,exclude=exclude.table,harmo_df=data_harm,harmo_exclude=harmo_exclude, ld_mat = ld_mat))
    
  }
  
}
  


  
  
  





