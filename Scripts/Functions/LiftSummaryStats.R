


###############################################################################
# This simple function will lift genomic coordinates (chromosome and position) from
# one build to another, namely from GRCh37 to GRCh38 or from GRCh38 to GRCh37.
#
# Note that GRCh37 is sometimes called hg19 and GRCh38 is sometimes called hg38.
#
# Arguments:
#    data_in: the input file containing 'chr' and 'pos'. The variable 'chr' should be formatted as 
#             'chrXX'. 'pos' should be numeric.
#    direction: the direction of the lift. Choices are '37_to_38' (default) and '38_to_37'
#
# The function will output a data frame with an additional column with the new coordinates, either 
#     'pos_ch37' or 'pos_ch38' depending on the lift selection
#
###############################################################################


LiftSummaryStats <- function(data_in, direction="37_to_38") {
  
  #Packages
  require(dplyr)
  require(liftOver)
  # Data management
  using_data<-data_in
  #Check duplicates
  using_data<- using_data %>% group_by(chr, pos) %>% dplyr::mutate(index = row_number())
  if (max(using_data$index>1)) {
    cat("There are duplicate combinations of chr and pos in the data; this should not affect the process")
  }
  using_data <- using_data %>% group_by(chr, pos) %>% dplyr::slice(1)
  # Fix value types
  using_data$chr <-as.character(using_data$chr)
  using_data$pos <-as.numeric(using_data$pos)
  #Check missing
  if (any(is.na(using_data$chr))|any(is.na(using_data$pos))) {
    cat("Some values for chr or pos are missing; Check your data. these will be ignored.")
    using_data<-using_data %>% filter(!is.na(pos)&chr!="")
  }
  df.tolift <- using_data %>% ungroup() %>%  transmute(chrom=chr, start=pos, end=pos,id = row_number())
  gr1 <- GRanges(df.tolift)
  

  if (direction=="37_to_38") {
    
    chain<- BiocIO::import("K:/isise/Procardis Topics/Proteomics QTLs/Data/ReceivedData/hg19ToHg38.over.chain")
    df.lifted <- as.data.frame(liftOver(gr1, chain))
    df.comb <- df.tolift %>% left_join(df.lifted,by="id") %>% 
      dplyr::select(chr=chrom, pos=start.x, pos_ch38 = start.y)
    data_out <- data_in   %>% left_join(df.comb, by=c("chr","pos"))
  }
  
  if (direction=="38_to_37") {
    
    chain<- BiocIO::import("K:/isise/Procardis Topics/Proteomics QTLs/Data/ReceivedData/hg38ToHg19.over.chain")
    df.lifted <- as.data.frame(liftOver(gr1, chain))
    df.comb <- df.tolift %>% left_join(df.lifted,by="id") %>% 
                            dplyr::select(chr=chrom, pos=start.x, pos_ch37 = start.y)
    data_out <- data_in   %>% left_join(df.comb, by=c("chr","pos"))
  }
  
  return(data_out)
  
  
  
}  







