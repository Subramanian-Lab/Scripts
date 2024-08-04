# ------------------------------------------------------------------------------
# Title: Calculating normalised signal for each chromosome position
# Author: Dr. Viji Subramanian
# Date: 21 June 2017
# Modified by: Devanarayanan P
# Modified on: 30 July 2024
# ------------------------------------------------------------------------------

generating_normalised_signal = function(genome_name = "SK1Yue", bedgraph_file_path, tw = 300) {
  
  message("Bedgraph file being imported ....")
  signal_data <- import_bedGraph(bedgraph_file_path)
  
  message("Retrieving genome coordinates for genome: ", genome_name)
  if (genome_name %in% c("SK1Yue", "S288CYue", "sacCer3", "SK1", "SK1_S288CYue")) {
    
    genome_info <- hwglabr2::get_chr_coordinates(genome_name)
  } else {
    
    if (is.null(gff_file_path)) {
      
      stop("gff file should be provided for the specified genome!!")
    }
    genome_info <- get_chromosome_cordinates(gff_file_path, genome_name)
  }
  
  # Calculating average signal and normalising signal data
  message("Computing average signal...")
  avrg <- function(x) (sum(GenomicRanges::width(x) * GenomicRanges::score(x), na.rm=TRUE) / sum(GenomicRanges::width(x)))
  genome_wide_mean <- avrg(signal_data)
  signal_data$score <- signal_data$score/genome_wide_mean
  
  # Sorting sequence levels in sample consistently and sorting genome info the same way
  sorted_signal <- sort(GenomeInfoDb::sortSeqlevels(signal_data))
  genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
  
  # Add chromosome length info to signal object
  GenomeInfoDb::seqlengths(sorted_signal) <- GenomeInfoDb::seqlengths(genome_info)
  
  # Compute 100-bp tiling windows
  bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(sorted_signal),
                                    tilewidth=tw, cut.last.tile.in.chrom=TRUE)
  
  # Get signal as "RleList"; the signal is stored in the "score" metadata column
  score <- GenomicRanges::coverage(sorted_signal, weight="score")
  
  # Compute average signal per tile
  bins <- GenomicRanges::binnedAverage(bins, score, "binned_score")
  
  # Get positions as the midpoints of the intervals
  positions <- bins@ranges@start + floor(bins@ranges@width / 2)
  
  # Make data frame (convert positions to Kb; signal is the binned score)
  message("Generating dataframe .... ")
  df <- data.frame(chr=bins@seqnames, position=positions / 1000, signal=bins$binned_score)
  df$chrSize <- GenomeInfoDb::seqlengths(bins)[df$chr]
  df <- dplyr::arrange(df, chrSize)
  
  return(df)
}