
rm(list=ls())

##################################################################################
# Install packages
#################################################################################

# CRAN

install.packages("rio")
install.packages("tidyverse")
install.packages("data.table")
install.packages("MendelianRandomization")
install.packages("testthat")
install.packages("grid")
install.packages("gridExtra")
install.packages("ieugwasr")
install.packages("patchwork")
install.packages("remotes")
install.packages("survey")
install.packages("markdown")
install.packages("pacman")

# From source

remotes::install_github("MRCIEU/TwoSampleMR")
install.packages("https://yanglab.westlake.edu.cn/software/gsmr/static/gsmr_1.1.1.tar.gz",repos=NULL,type="source")


##################################################################################
# Load packages
#################################################################################

if (!require(pacman)) install.packages("pacman")
p_load(rio, tidyverse,data.table, MendelianRandomization,
testthat, grid,gridExtra,ieugwasr,patchwork,TwoSampleMR,gsmr )



##################################################################################
# Load function and dependencies
#################################################################################

# Load main function

devtools::source_url("https://raw.githubusercontent.com/adamkvonende/MRpipeline/master/Scripts/perform_mr.R")
#source("K:/isise/Procardis Topics/Proteomics QTLs/Analyses/MRpipeline/Scripts/perform_mr.R")

# Load supporting functions


devtools::source_url("https://raw.githubusercontent.com/adamkvonende/MRpipeline/master/Scripts/Functions/LiftSummaryStats.R")
devtools::source_url("https://raw.githubusercontent.com/adamkvonende/MRpipeline/master/Scripts/Functions/Map_rsid_from_chrpos.R")


source("K:/isise/Procardis Topics/Proteomics QTLs/Analyses/MRpipeline/Scripts/perform_mr.R")

# Load supporting functions

source("K:/isise/Procardis Topics/Proteomics QTLs/Analyses/MRpipeline/Scripts/Functions/LiftSummaryStats.R")
source("K:/isise/Procardis Topics/Proteomics QTLs/Analyses/MRpipeline/Scripts/Functions/Map_rsid_from_chrpos.R")


###############################################################################
# Load exposure and outcome data
###############################################################################

# Load exposure data: summary statistics for Apolipoprotein C3 (ApoC3)

setwd("K:/isise/Procardis Topics/Proteomics QTLs/Analyses/MRpipeline/Data")


ss_apoc3 <- rio::import("https://zenodo.org/record/8090662/files/apoc3_summary_stats_pless001.csv?download=1")

# Load outcome data: summary statistics for coronary heart disease (CHD)

ss_chd <- rio::import("https://zenodo.org/record/8090662/files/chd_summary_stats.csv?download=1") 
ss_apoc3 <- rio::import("apoc3_summary_stats_pless001.csv")

# Load outcome data: summary statistics for coronary heart disease (CHD)

ss_chd <- rio::import("chd_summary_stats.csv") 


###############################################################################
# Run an MR analysis of ApoC3 and CHD using all variants significantly (p < 5x10^-8) associated with 
# ApoC3 concentration (all methods)
###############################################################################

# Note: we run with fewer PRESSO bootstraps to save time.



results1 <- perform_mr(ss_exposure=ss_apoc3, ss_outcome=ss_chd, cis = F, 
                          maf_thresh=0.01, p_thresh=5*10^-8, 
                          harmo_action=3, # drop palindromic
                          clump=T,clump_kb=250, clump_r2 = 0.01,clump_pop="EUR", 
                          name_exp="Apoc3", name_outcome="CHD",
                          which_methods="all", 
                          ld_mat = "try_gen", ld_mat_pop="EUR",
                          presso_nboot = 1000,which_plots="all", 
                          outcome_binary=T)
       

# View outputs


results1$results %>% View()
results1$het %>% View()
results1$exclude %>% View()
results1$harmo_exclude %>% View()



results1$results %>% View()
results1$het %>% View()
results1$exclude %>% View()
results1$for_plot %>% View()
results1$for_plot

results.all$results %>% View()
results.all$het %>% View()
results.all$exclude %>% View()
results.all$harmo_exclude %>% View(



###############################################################################
# Run a cis-MR analysis of ApoC3 and CHD using all variants significantly (p < 5x10^-8) associated 
# with ApoC3 concentration that reside within 100kb of the APOC3 gene 
# (basic methods only)
###############################################################################



results2<- perform_mr(ss_exposure=ss_apoc3, ss_outcome=ss_chd, cis = T, 
                          gene_chr=11, gene_start=116700422,gene_end=116703788,
                          flank_kb=100,maf_thresh=0.01, p_thresh=5*10^-8, 
                          harmo_action=3, clump=F,clump_kb=250, clump_r2 = 0.001,clump_pop="EUR", 
                          name_exp="Apoc3", name_outcome="CHD",
                          which_methods="basic", which_plots="all", 
                          outcome_binary=T)





results2$results %>% View()
results2$het %>% View()
results2$harmo_exclude %>% View()


####################
# Perform clumping only.
# This may be useful if you want to create an independent set of variants irrespective
# of a specific outcome
###############################################################################

clumped_list <- perform_mr(ss_exposure=ss_apoc3, 
                            p_thresh=5*10^-8, 
                           clump_only=T,clump_kb=250, clump_r2 = 0.1,clump_pop="EUR")
                                           
