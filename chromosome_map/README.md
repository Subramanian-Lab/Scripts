# Chromosome map

The code provided can be used to plot chromosome maps for *S. Cerevisiae* like shown in [Subramanian VV. et al. 2019](https://www.nature.com/articles/s41467-019-08875-x) for custom genome. 

## Included files

### Functions: 
- `generating_normalised_signals.R` -> Calculating normalised signal for each chromosome position.  
- `get_chromosome_cordinates` -> Creating a GRanges object containing centromere data from a custom .gff file 
- `generating_centromere_data.R` -> Calculating centromere midpoints for the chromosomal plot. 

### Examples case:
- `plotting_chromosome_map.R` -> shows an example code for using the functions to plot the chromosome map


## Inputs
1. Name of the strain:
    - If the strain name is "SK1Yue", "S288CYue", "sacCer3", "SK1", "SK1_S288CYue", then the chromosome centromere data are obtained from the `get_chr_cordinates` function from the [hwglabr2](https://github.com/hochwagenlab/hwglabr2) R-package. 
2. If the strain is not listed above, a custom .gff file of the genome can be used to generate the centromere data using `get_chromosome_cordinates.R` file.
3. The bedgraph file containing the positions of sequences and its scores to plot.

## Output
The output will be a plot which can be exported out.   
