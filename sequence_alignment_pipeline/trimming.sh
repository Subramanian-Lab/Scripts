#!/bin/zsh

## This code trims and remove the adapter sequences from the raw reads and removes low quality bases

raw_file_dir="../raw_files"

for file in "$raw_file_dir"/dataset*; do
    echo $file
    file_base=$(basename "$file" .fastq)
    trim_galore -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --length 15 -q 33 $file 
    echo "Trimmed reads in $file_base : $count"
done
