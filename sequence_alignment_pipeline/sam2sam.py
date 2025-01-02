"""
This code reassigns alignment score by penalising soft masking ends in the CIGAR strings, and then segregate into
unique and multi mapped read files by making use of NH tag of SAM file
"""

# Importing necessary libraries
import sys
import pysam
from tqdm import tqdm

# Clipping values in pysam
bam_soft_clip = 4
bam_hard_clip = 5
bam_match = 0

# creating a cigar dictionary to create cigar strings
cigar_dictionary = {bam_soft_clip: "S",bam_match: "M", bam_hard_clip: "H"}

# Function to convert CIGAR tuple into a CIGAR string (from pysam output)
def formatCigar(cigar):
    cigarStr=""
    for (code,num) in cigar:
        try:
            cigarStr+=str(num)+cigar_dictionary[code]
        except:
            print(cigar, code, num)
            raise
    return cigarStr

# Function to get the arguments given to the gmapper-ls in sequence alignment
def getGmapperOpts(cmd):
    F=cmd.split()
    parse=[(i,x.replace("-","")) for i,x in enumerate(F) 
           if x[0]=="-" and x[1] not in "-0123456789"]
    ret=dict([(x,F[i+1]) for i,x in parse]) 
    ret['h']=int(ret['h'])
    ret['i']=int(ret['i'])
    ret['m']=int(ret['m'])
    return ret

# Obtaining SAM alignment positions and strands
def fmtSAMalignPos(aa):
    if aa.flag & 16:
        return [aa.aend,aa.pos+1,"-"]
    else:
        return [aa.pos+1,aa.aend,"+"]

# Initialising input and output files
sam = pysam.AlignmentFile(sys.argv[1])
opts = getGmapperOpts(sam.header["PG"][0]["CL"])
unique_output_file = pysam.AlignmentFile(f"{sys.argv[1][-20:-4]}_unique.sam", "w", header=sam.header)
multi_output_file = pysam.AlignmentFile(f"{sys.argv[1][-20:-4]}_multi.sam", "w", header=sam.header)
total_file = pysam.AlignmentFile(f"./{sys.argv[1][-20:-4]}_total.sam", "w", header=sam.header)

# Counters to count for unique and multi mapped reads
counter_unique = 0
counter_multi = 0
total_no = 0

print(f"Processing {sys.argv[1]}")

# Looping through each reads to remove soft clipped ends
for si in tqdm(sam):
    if si.is_unmapped:
        continue
    total_no += 1

    left_clip = 0 if si.cigar[0][0] != bam_soft_clip else si.cigar[0][1]
    left_clip_seq = si.seq[:left_clip]
    left_clip_mismatch = sum([x.upper() != "C" for x in left_clip_seq])

    right_clip = 0 if si.cigar[-1][0] != bam_soft_clip else si.cigar[-1][1]
    right_clip_seq = si.seq[-right_clip:]
    right_clip_mismatch = sum([x.upper() != "G" for x in right_clip_seq])

    old_alignment_score = si.opt("AS")
    new_alignment_score = opts["i"] * (left_clip_mismatch + right_clip_mismatch) + old_alignment_score

# Segregating reads to unique and multi based on NH value
    if new_alignment_score >= opts["h"]:
        si.set_tag("AS", new_alignment_score)
        total_file.write(si)
        if si.opt("NH") == 1:   
            counter_unique += 1
            unique_output_file.write(si)
        else:
            counter_multi += 1
            multi_output_file.write(si)

# Printing results
print(f"Unique reads: {counter_unique}")
print(f"Multi: {counter_multi}")
print(f" Total: {total_no}")
