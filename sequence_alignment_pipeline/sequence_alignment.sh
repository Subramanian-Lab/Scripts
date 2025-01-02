#!/bin/bash

# This code aligns the split mixed sequences with the reference genome with the help of SHRiMP sequence aligner

SHRIMP=./SHRiMP_2_1_1/bin/gmapper-ls
GENOME="SK1.genome.fa"
QUERY_FILE="$1"
BASE=$(basename "$QUERY_FILE" .fq)

# Function for testing file existence
test_file_existence() {
    file_name=$1
    if [ ! -f $file_name ]
    then
        echo "The file $file_name does not exist"
    else
        echo "Test pass"
    fi
}

# Testing genome file and query file
test_file_existence $GENOME
test_file_existence $QUERY_FILE

# Sequence alignment
$SHRIMP -N 15 -U -g -1000 -q -1000 \
    -m 10 -i -20 -h 100 -r 50% \
    -n 1 -s 1111111111,11110111101111,1111011100100001111,1111000011001101111 \
   -o 1001 -Q -E --sam-unaligned --strata \
  $QUERY_FILE $GENOME > "${BASE}.sam" 2>"${BASE}.log" 

echo "$QUERY_FILE aligned and results saved... "
