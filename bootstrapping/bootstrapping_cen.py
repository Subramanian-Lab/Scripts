# Importing necessary packages
from collections import defaultdict
import random
import matplotlib.pyplot as plt
from datetime import datetime
import time, sys
from tqdm import tqdm

# File paths
bedgraph_file_path = "../../../mapping_project/zip3_del/peakcalling/dataset1_bdgcmp.bdg"
gff_file_path = "../../../mapping_project/zip3_del/peakcalling/SK1.all_feature.gff"

# Reading the bedgraph file
with open(bedgraph_file_path, "r") as f:
    bedgraph_list = [line.strip().split('\t') for line in f]

bedgraph_dictionary = defaultdict(list)

for i in bedgraph_list:
    if len(i) != 4:
        print("Rows not complete")
        exit()
    else:
        bedgraph_dictionary[i[0]].append(list(map(float, i[1:])))

# Reading the bedgraph dictionary into a wiggle dictionary with chromosome keys
wiggle_dictionary = defaultdict(list)

for chr in bedgraph_dictionary.keys():
    for row in range(len(bedgraph_dictionary[chr])):
        if row != 0:
            # Making the interval continuous
            if bedgraph_dictionary[chr][row][0] < bedgraph_dictionary[chr][row - 1][1]:
                bedgraph_dictionary[chr][row][0] = bedgraph_dictionary[chr][row - 1][1]

            # Writing per base signal to wiggle dictionary
            if bedgraph_dictionary[chr][row][2] != 0:
                positions = list(range(int(bedgraph_dictionary[chr][row][0]) + 1, int(bedgraph_dictionary[chr][row][1]) + 1))
                wiggle_dictionary[chr] += [[position, bedgraph_dictionary[chr][row][2]] for position in positions]

# Calculating genome average signal from the wiggle dictionary
data_genome = defaultdict(list)
for chr in wiggle_dictionary.keys():
    signal = [i[1] for i in wiggle_dictionary[chr]]
    data_genome[chr].append([sum(signal), len(signal)])

genome_signal = sum(i[0][0] for i in data_genome.values())
genome_length = sum(i[0][1] for i in data_genome.values())

genome_average_signal = genome_signal / genome_length

# Converting bed or gff file to universal format
with open(gff_file_path, "r") as f:
    input = [row.strip().split("\t") for row in f if not row.startswith("#")]

# Taking only the centromere data
input = [row for row in input if "centromere" in row]

if len(input) == 9:
    for i in range(len(input)):
        chr = input[i][0]
        # Reformating to [chr, start, end, signal, id]
        input[i] = [chr, int(input[i][3]), int(input[i][4]), input[i][8], input[i][6]]

# Creating signal average for centromere
data_real = defaultdict(list)
extend = 2000

for cen in input:
    chr = cen[0]
    mid = round((cen[1] + cen[2]) / 2)
    start = mid - (extend/2)
    end = mid + (extend/2)
    cen_signal_list = [i[1] for i in wiggle_dictionary[chr] if i[0] >= start and i[0] <= end]
    data_real[chr].append([sum(cen_signal_list), len(cen_signal_list)])

cen_signal = sum([i[0][0] for i in data_real.values()])
cen_length = sum([i[0][1] for i in data_real.values()])

# Normalising with genome average signal
cen_ratio = (cen_signal/cen_length) / genome_average_signal

chr_order = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"]

# Writing the cen signal result to file
with open("outputfile.txt", "w") as f:
    f.write("\t".join(['region', 'total_signal', 'number_of_bases']) + "\n")
    for chr in chr_order:
        f.write(f"{chr}\t{str(data_real[chr][0][0])}\t{str(data_real[chr][0][1])}\n")
    f.write(f"genome\t{str(genome_signal)}\t{str(genome_length)}\n")
    f.write(f"ratio:\t{str(cen_ratio)}")

# Random sampling signals
num_realisation = 5000
data_rand = defaultdict(list)
rand_signal = []
rand_length = []

for i in tqdm(range(num_realisation)):
    for chr in wiggle_dictionary.keys():
        random_position = random.randint(len(wiggle_dictionary[chr]) - 1)
        start = wiggle_dictionary[chr][random_position][0]
        end = start + extend

        if wiggle_dictionary[chr][-1][0] < end:
            new_end = end - wiggle_dictionary[chr][-1][0] + 1
            signal_list = [j[1] for j in wiggle_dictionary[chr] if j[0] >= start or j[0] <= end]
            data_rand[chr].append([sum(signal_list), len(signal_list)])
        else:
            signal_list = [j[1] for j in wiggle_dictionary[chr] if j[0] >= start and j[0] <= end]
            data_rand[chr].append([sum(signal_list), len(signal_list)])
            
    rand_signal.append(sum([data_rand[key][i][0] for key in data_rand.keys()]))
    rand_length.append(sum([data_rand[key][i][1] for key in data_rand.keys()]))

rand_ratio = [(rand_signal[i]/rand_length[i]) / genome_average_signal for i in range(len(rand_signal))]

with open("rand_ratio.txt", "w") as f:
    f.write("\n".join(map(str, rand_ratio)))

plt.hist(rand_ratio, bins=20)
plt.axhline(cen_ratio, color="b", linestyle="dashed", linewidth=2)
plt.title("Bootstrap of ratios")
plt.xlabel("Feature average / average genome signal")
plt.ylabel("Count")
plt.show()
