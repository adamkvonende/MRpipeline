


###############################################################################
# This simple function will attempt to find rsids based on chromsosome and position for different
# builds (GRCh37 or GRCh38)
# Note that GRCh37 is sometimes called hg19 and GRCh38 is sometimes called hg38.
#
# Arguments:
#    data_in: the input file containing 'chr' and 'pos'. The variable 'chr' should be formatted as 
#             'XX', e.g. 1, 2,3, 15, etc
#    build: the build on which the chrosome and pos coordinates are based ("ch37" or "ch38")
#
# The function will output a data frame with an additional column with the mapped rsid
#
# 16 March 2023
#
# ch37 is based on SNPlocs.Hsapiens.dbSNP144.GRCh37
# ch38 is based on SNPlocs.Hsapiens.dbSNP155.GRCh38
#
# You will need to install the following packages in your local directory. 
#  BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37",lib="C:/Users/adamv/AppData/Local/R/win-library/4.2")
#  BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38",lib="C:/Users/adamv/AppData/Local/R/win-library/4.2")
#  Use .libPaths() to identify your local directory
###############################################################################


# Load packages
# library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
# library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
# library(dplyr)



Map_rsid_from_chrpos <- function(data_in, build="ch37") {
  
  start_time <- Sys.time()
  # Data management
  # Fix value types
  data_in$chr <-as.character(data_in$chr)
  data_in$pos <-as.numeric(data_in$pos)
  using_data<-data_in
  #Check duplicates
  using_data<- using_data %>% group_by(chr, pos) %>% dplyr::mutate(index = row_number())
  if (max(using_data$index>1)) {
    cat("There are duplicate combinations of chr and pos in the data; this should not affect the process")
  }
  using_data <- using_data %>% group_by(chr, pos) %>% dplyr::slice(1)
 
  #Check missing
  if (any(is.na(using_data$chr))|any(is.na(using_data$pos))) {
    cat("Some values for chr or pos are missing; These will not be mapped\n\n")
    using_data<-using_data %>% filter(!is.na(pos)&chr!="")
  }
  
  using_data <-using_data %>% transmute(chr, start=pos, end=pos)
  map.ranges <- makeGRangesFromDataFrame(using_data)
  
  if (build == "ch37") {
    out <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, map.ranges )
  } 
  if (build == "ch38"){
    out <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, map.ranges )
  }
  out.df <- as.data.frame(out) %>% transmute(chr = seqnames,pos, rsid_mapped=RefSNP_id)
  cat(nrow(out.df), "rsids were successfuly mapped\n\n")
  
  # merge
  data_out <- data_in %>% left_join(out.df, by=c("chr","pos"),multiple="all")
  end_time <- Sys.time()
  time_taken = end_time - start_time
  mins <- as.integer(time_taken / 60)
  cat("rsid mapping took", mins, "minutes\n\n")
  return(data_out)
}
  
#out.mapped <- Map_rsid_from_chrpos(test3,build="ch38")
  