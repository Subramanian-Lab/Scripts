"""
This code removes the 5 bases from the start of each read (diffusion primers from the library preparation)
and adds it into the query name for downstream identification
"""
# Importing necessary libraries
import sys
import Bio.SeqIO
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Specifying input and output directories
input_dir = "../trimmed_files"
output_dir = "../splitmixed"

# Looping through each read to modify the sequences and writing into the output file
for file in os.listdir(input_dir):
    if file.endswith(".fq"):
        read_file = os.path.join(input_dir, file)
        print(read_file)
        modified_sequences = []
        for seq in Bio.SeqIO.parse(read_file, "fastq"):
            newSeq=seq[5:]
            newSeq.id = seq.id + ":" + str(seq.seq[:5])
            newSeq.name = ""
            newSeq.description = ""
            newSeq.seq = Seq(str(seq.seq)[5:])
            modified_sequences.append(newSeq)
        output_file = os.path.join(output_dir, file)
        print(output_file)
        with open(output_file, "w") as out_handle:
            Bio.SeqIO.write(modified_sequences, out_handle, "fastq")
