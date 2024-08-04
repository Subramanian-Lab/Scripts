# ------------------------------------------------------------------------------
# Title: Calculating centromere midpoint for chromosome plot
# Author: Dr. Viji Subramanian
# Date: 21 June 2017
# Modified by: Devanarayanan P
# Modified on: 30 July 2024
# ------------------------------------------------------------------------------

generating_centromere_data = function(genome_name = "SK1Yue", gff_file_path = NULL) {
  if (genome_name %in% c("SK1Yue", "S288CYue", "sacCer3", "SK1", "SK1_S288CYue")) {
    centromeres <- hwglabr2::get_chr_coordinates(genome_name)
  } else {
    if (is.null(gff_file_path)) {
      stop("gff file should be provided for the specified genome...")
    }
    centromeres <- get_chromosome_cordinates(gff_file_path, genome_name)
  }
  
  positions <- GenomicRanges::start(centromeres) + floor(GenomicRanges::width(centromeres) / 2)
  positions
  df2 <- data.frame(chr=centromeres@seqnames, position=positions / 1000)
  df2
  df2$chrSize <- GenomeInfoDb::seqlengths(centromeres)[df2$chr]
  df2 <- dplyr::arrange(df2, df2$chrSize)
  
  return(df2)
}