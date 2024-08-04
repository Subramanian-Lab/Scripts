# Sourcing files
source("get_chromosome_cordinates.R")
source("generating_centromere_data.R")
source("generating_normalised_signal.R")

# Importing necessary libraries
library(GenomeInfoDb)
library(GenomicRanges)
library(tidyverse)
library(hwglabr2)

# file paths for bedgraph and gff
gff_file_path <- file.path("file.gff")
bedgraph_file_path <- file.path("file.bedgraph")

# Provide the strain name
strain_name <- "StrainName" 

# Obtaiing Centromere data
df <- generating_centromere_data(strain_name, gff_file_path)

# Obtaining signal data
signal_data <- generating_normalised_signal(strain_name, bedgraph_file_path)

# Plotting the signal and centromeres
p1 <- ggplot(signal_data, aes(x=position, y=signal)) + 
  geom_line(colour="blue") + 
  geom_point(data =df, aes(x=position, y=-1.0), color='black', size=1.9, pch=17)

# Facet gridding multiple chromosomes
p1 + facet_grid(chrSize ~ .) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text.y = element_blank(), panel.background = element_blank(), panel.spacing = unit(0, "lines"))

