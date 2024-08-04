# ------------------------------------------------------------------------------
# Title: Creating a GRanges object containing centromere data for a custom
#       genome.
# Author: Devanarayanan P
# Date: 30 July 2024
# Based on: hwlabgr2 Github page (data_internal.R)
# (https://github.com/hochwagenlab/hwglabr2/blob/master/data-raw/data_internal.R)
# ------------------------------------------------------------------------------

# Function to add the name of the genome as attribute to the GRanges
add_genome_name_to_GR <- function(g_ranges_object, name = 'SK1Yue') {
  if (is.null(seqinfo(g_ranges_object))) {
    stop("The GRanges object does not have sequence information.")
  }
  
  # Set the genome name
  GenomeInfoDb::genome(seqinfo(g_ranges_object)) <- name
  
  return(g_ranges_object)
}

# Function for creating G Ranges object for a custom .gff file (centromere pos)

get_chromosome_cordinates = function(gff_file_path, genome_name) {
  
  # Obtaining the centromere data from .gff file
  genome_gff <- rtracklayer::import.gff(gff_file_path)
  genome_cen <- genome_gff[genome_gff$type == "centromere"]
  
  # Dropping the chrMIto and 2-micron data to isolate the chromosome data
  genome_cen <- GenomeInfoDb::dropSeqlevels(genome_cen, c("chrMito", "2-micron"))
  GenomicRanges::mcols(genome_cen) <- NULL
  
  # Adding the chromosome lengths as an attribute to the GRanges
  chr_len <- genome_gff[genome_gff$type == "chromosome"][1:16]
  chr_names <- as.character(seqnames(chr_len)) 
  chr_len <- setNames(width(chr_len), chr_names)
  GenomeInfoDb::seqlengths(genome_cen) <- chr_len

  # Adding the genome name as an attribute to the GRanges object
  genome_cen <- add_genome_name_to_GR(genome_cen, name = genome_name)

  return(genome_cen)
}