# Import necessary packages
from collections import defaultdict
import random
from tqdm import tqdm

# File paths
bedgraph_file_path = "../../../mapping_project/zip3_del/peakcalling/dataset1_bdgcmp.bdg"
gff_file_path = "../../../mapping_project/zip3_del/peakcalling/SK1.all_feature.gff"

# Reading bedgraph file
def read_bedgraph(input_file):
    with open(input_file, "r") as f:
        bedgraph_list = [line.strip().split('\t') for line in f]

    # Creating a dictionary with chromosome order
    bedgraph_dictionary = defaultdict(list)

    for row in bedgraph_list:
        if len(row) != 4:
            print("Rows not complete")
            exit()
        else:
            bedgraph_dictionary[row[0]].append(list(map(float, row[1:])))

    # Sorting chromosome numerically
    for key in bedgraph_dictionary.keys():
        bedgraph_dictionary[key].sort()
    
    # Converting it into a wiggle dictionary
    wiggle_dictionary = defaultdict(list)
    
    for chr in bedgraph_dictionary.keys():
        for row in range(len(bedgraph_dictionary[chr])):
            if row != 0:
                # Making intervals continous
                if bedgraph_dictionary[chr][row][0] < bedgraph_dictionary[chr][row - 1][1]:
                    bedgraph_dictionary[chr][row][0] = bedgraph_dictionary[chr][row - 1][1]
                
                # assigning signals to invididual bases
                if bedgraph_dictionary[chr][row][2] != 0:
                    positions = list(range(int(bedgraph_dictionary[chr][row][0]) + 1, int(bedgraph_dictionary[chr][row][1]) + 1))
                    wiggle_dictionary[chr] += [[position, bedgraph_dictionary[chr][row][2]] for position in positions]
                
    return wiggle_dictionary

# Reading BED / GFF file with feature
def read_bed_gff_file(feature_file):
    with open(feature_file, "r") as f:
        feature_input = [row.strip().split('\t') for row in f if not row.startswith("#")]
    
    # to extract only lines containing centromere data
    feature_input = [row for row in feature_input if "centromere" in row]

    # For GFF file
    if len(feature_input[0]) == 9:
        for i in range(len(feature_input)):
            chr = feature_input[i][0]

            # Ordering as [chr, start, end, strand, id]
            feature_input[i] = [chr, int(feature_input[i][3]), int(feature_input[i][4]), feature_input[i][8], feature_input[i][6]]

    # For simple BED file
    elif len(feature_input[0]) == 3:
        for i in range(len(feature_input)):
            chr = feature_input[i][0]

            feature_input[i] = [chr, int(feature_input[i][1]), int(feature_input[i][2], f"id{str(i)}", "+" )]
        
    return feature_input

# Generating ratio from the feature data
def real_data_ratio(bed_file, wiggle_dictionary, extend):
    extend = int(extend)
    centromere_data = read_bed_gff_file(bed_file)

    data_real = defaultdict(list)

    # Averaging signal around an extended region from the midpoint
    for cen in centromere_data:
        chr = cen[0]
        mid = round((cen[1] + cen[2]) / 2)
        start = mid - (extend / 2)
        end = mid + (extend / 2)
        cen_signal_list = [row[1] for row in wiggle_dictionary[chr] if row[0] >= start and row[0] <= end]
        data_real[chr].append([sum(cen_signal_list), len(cen_signal_list)])
    
    total_cen_signal = sum([i[0][0] for i in data_real.values()])
    total_cen_length = sum([i[0][1] for i in data_real.values()])

    return data_real, total_cen_length, total_cen_signal

# Taking ratio from random sampling
def random_data_ratio(num_realisation, wiggle_dictionary, extend):
    data_random = defaultdict(list)
    total_random_signal = []
    total_random_length = []

    for i in tqdm(range(int(num_realisation))):
        for chr in wiggle_dictionary.keys():
            random_index = random.randint(0, len(wiggle_dictionary[chr]) - 1)
            start = wiggle_dictionary[chr][random_index][0]
            end = start + extend

            # wrapping around if the end point is greater than genome length
            if wiggle_dictionary[chr][-1][0] < end:
                new_end = end - wiggle_dictionary[chr][-1][0] + 1
                signal_list = [j[1] for j in wiggle_dictionary[chr] if j[0] >= start or j[0] <= end]
                data_random[chr].append([sum(signal_list), len(signal_list)])

            else:
                signal_list = [j[1] for j in wiggle_dictionary[chr] if j[0] >= start and j[0] <= end]
                data_random[chr].append([sum(signal_list), len(signal_list)])

        total_random_signal.append(sum([data_random[key][i][0] for key in data_random.keys()]))
        total_random_length.append(sum([data_random[key][i][1] for key in data_random.keys()]))

    return total_random_signal, total_random_length

# Taking genome average signal for normalisation
def calculate_genome_average_signal(wiggle_dictionary):
    data_genome = defaultdict(list)

    for chr in wiggle_dictionary.keys():
        signal_list = [i[1] for i in wiggle_dictionary[chr]]
        data_genome[chr].append([sum(signal_list), len(signal_list)])
    
    total_genome_signal = sum([chr_signal[0][0] for chr_signal in data_genome.values()])
    total_genome_length = sum([chr_signal[0][1] for chr_signal in data_genome.values()])

    return total_genome_length, total_genome_signal

def main(bedgraph_file, extend, num_realisation, feature_file, output_tag):
    wiggle_dictionary = read_bedgraph(bedgraph_file)

    # Calculating genome average signal
    total_genome_signal, total_genome_length = calculate_genome_average_signal(wiggle_dictionary)
    genome_average_signal = total_genome_signal / total_genome_length

    # Calculating average signal at centromere (normalised)
    data_real, total_cen_signal, total_cen_length = real_data_ratio(feature_file, wiggle_dictionary, extend)
    cen_ratio = (total_cen_signal / total_cen_length) / genome_average_signal

    # Writing results to output file
    chr_order = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"]

    with open(f"{output_tag}_cendata_ext-{str(extend)}.txt", "w") as f:
        f.write("\t".join(["region", "total_signal", "num_bases"]) + "\n")
        for chr in chr_order:
            f.write(f"{chr}\t{str(data_real[chr][0][0])}\t{str(data_real[chr][0][1])}\n")
        f.write(f"genome\t{str(total_genome_signal)}\t{str(total_genome_length)}\n")
        f.write(f"ratio:\t{str(cen_ratio)}")

    # Calculating average signal by random sampling
    total_random_signal, total_random_length = random_data_ratio(num_realisation, wiggle_dictionary, extend)
    rand_ratio = [(total_random_signal[i] / total_random_length[i]) / genome_average_signal for i in range(len(total_random_signal))]

    # Writing results to output file
    with open(f"{output_tag}_randratio_num-real-{str(num_realisation)}_ext-{str(extend)}.txt", "w") as f:
        f.write("\n".join(map(str, rand_ratio)))

main(bedgraph_file_path, 2000, 5000, gff_file_path, "zip3_del")
