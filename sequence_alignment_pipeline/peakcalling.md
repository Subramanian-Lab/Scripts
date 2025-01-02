# Instructions to do peak calling using Picard Tools and MACS3

## <a name="contents"></a> Contents
1. [List of files](#filelist)
2. [Order of Excecution](#order)
3. [Code Credits](#credits)

## <a name="filelist"></a> List of Files
1. trimming.sh -> Trims of adapters and low quality bases.
2. splitMixer.sh -> Trims of 5 bases from the start of each read and adds it to the query name for downstream identification.
3. sequence_alignment.sh -> splitmixed reads are sequence aligned.
4. alignment.sh -> runs the sequence_alignment.sh scripts for each of the trimmed reads.
5. sam2sam.py -> penalising soft mapped reads and reassigning mapping scores. Also seperates unique and multimapped reads.

## <a name="order"></a> Order of Execution
1. trimming.sh
2. splitMixer.sh
3. alignment.sh
4. sam2sam.py

## <a name="credits"></a> Code Credits
[D. Thacker et al. 2011](https://pubmed.ncbi.nlm.nih.gov/24717437)
