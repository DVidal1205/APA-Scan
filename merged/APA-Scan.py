import time
import pandas as pd
import methods
import glob
import os
import configparser
from multiprocessing import cpu_count

# Set the flags for the species, input, output, and PASseq directories
(speciesFlag, inputFlag, outFlag, pasFlag) = (0, 0, 0, 0)

# Begin the benchmark timer
startTime = time.time()


def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, "*")))]


# Read in the configuration values
config = configparser.ConfigParser()
config.read("configuration.ini")
input1_dir = config["INPUT_RNAseq"]["input1"]
input2_dir = config["INPUT_RNAseq"]["input2"]
if input1_dir[-1] == "/":
    input1_dir = input1_dir[:-1]
if input2_dir[-1] == "/":
    input2_dir = input2_dir[:-1]
pasSeq_dir1 = config["INPUT_PASseq"]["pas1"]
pasSeq_dir2 = config["INPUT_PASseq"]["pas2"]
if pasSeq_dir1[:-1] == "/":
    pasSeq_dir1 = pasSeq_dir1[:-1]
if pasSeq_dir2[:-1] == "/":
    pasSeq_dir2 = pasSeq_dir2[:-1]
output_dir = config["OUTPUT_FOLDER"]["output_dir"]
if output_dir[-1] != "/":
    output_dir += "/"
extended = config["Extended_3UTR"]["extended"]
all_events = config["All_events"]["All"]

############################################################
# NEW CODE: Handle cores input from the configuration file #
############################################################
parallel = True

cores = config["CORES"]["cores"]

# Cases where we run in sequential mode
if cores == "NULL":
    parallel = False
if cores == "1":
    print("Running with one core is slower than simply running the code in sequential mode. Running APA-Scan sequentially .")
    parallel = False

# Check for MAX case
if cores == "MAX":
    cores = cpu_count()

# Let the user know how many cores we are using and if we are running in parallel
if parallel:
    print("Running APA-Scan in parallel mode with", cores, "cores")
else:
    print("Running APA-Scan in sequential mode")

# Inform the user of the configuration settings
os.makedirs(output_dir, exist_ok=True)
inp_annotation = config["ANNOTATION"]["annotation"]
ref_genome = config["ANNOTATION"]["genome"]
g1_name, g2_name = os.path.basename(input1_dir), os.path.basename(input2_dir)
print("RNA-seq input 1 dir:", input1_dir)
print("RNA-seq input 2 dir:", input2_dir)
print("3'-end-seq input 1 dir:", pasSeq_dir1)
print("3'-end-seq input 2 dir:", pasSeq_dir2)
print("Output Dir:", output_dir)
print("Annotation:", inp_annotation, ref_genome, "\n")
print("Loading annotation information ...")
ann_df = pd.read_csv(inp_annotation, delimiter="\t")
chr_list = ann_df["chrom"].str.split("_", n=1, expand=True)
chromosomes = list(set(chr_list[0]) - set(["chrM", "chrUn"]))

# Clean the dataframe
if "#name" in ann_df.columns:
    ann_df.rename(columns={"#name": "name"}, inplace=True)

# Begin generating the read coverage files
print("Generating read coverage files from the RNA-seq data...")
for sample1 in os.listdir(input1_dir):
    ####################################################################
    # NEW CODE: Handle parallel mode for SamtoText if parallel is True #
    ####################################################################

    if sample1.endswith(".bam"):
        if parallel:
            methods.SamtoTextParallel(input1_dir, sample1, chromosomes, cores)
        else:
            methods.SamtoTextSequential(input1_dir, sample1, chromosomes)

for sample2 in os.listdir(input2_dir):
    ####################################################################
    # NEW CODE: Handle parallel mode for SamtoText if parallel is True #
    ####################################################################

    if sample2.endswith(".bam"):
        if parallel:
            methods.SamtoTextParallel(input2_dir, sample2, chromosomes, cores)
        else:
            methods.SamtoTextSequential(input2_dir, sample2, chromosomes)

# Specify the result filename
result_filename = "APA_Scan_" + g1_name + "_Vs_" + g2_name

# Process the RNA-seq data since we don't have 3'-end-seq data
if pasSeq_dir1 == "NULL" or pasSeq_dir2 == "NULL":
    s1_namelist = list_dirs(input1_dir)
    s2_namelist = list_dirs(input2_dir)

    print("Preparing result using RNA-seq data only")
    methods.Get_Signal_Positions(chromosomes, ann_df, ref_genome, output_dir, extended)
    if all_events == "all":
        methods.with_PAS_signal_all(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1_name, g2_name, output_dir, result_filename)
    else:
        methods.with_PAS_signal(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1_name, g2_name, output_dir, result_filename)

else:
    # Begin handling the 3-seq data
    print("Generating read coverage files for the 3'-end-seq data...")
    for sample1 in os.listdir(pasSeq_dir1):
        ####################################################################
        # NEW CODE: Handle parallel mode for SamtoText if parallel is True #
        ####################################################################

        if sample1.endswith(".bam"):
            if parallel:
                methods.SamtoTextParallel(pasSeq_dir1, sample1, chromosomes, cores)
            else:
                methods.SamtoTextSequential(pasSeq_dir1, sample1, chromosomes)

    for sample2 in os.listdir(pasSeq_dir2):
        ####################################################################
        # NEW CODE: Handle parallel mode for SamtoText if parallel is True #
        ####################################################################

        if sample2.endswith(".bam"):
            if parallel:
                methods.SamtoTextParallel(pasSeq_dir2, sample2, chromosomes, cores)
            else:
                methods.SamtoTextSequential(pasSeq_dir2, sample2, chromosomes)

    p1_namelist = list_dirs(pasSeq_dir1)
    p2_namelist = list_dirs(pasSeq_dir2)
    p1_name = pasSeq_dir1.split("/")[-1]
    p2_name = pasSeq_dir2.split("/")[-1]

    methods.Get_Peak_Positions(chromosomes, ann_df, pasSeq_dir1, pasSeq_dir2, p1_name, p2_name, output_dir, extended)
    if all_events == "all":
        methods.with_PA_peaks_all(chromosomes, input1_dir, input2_dir, g1_name, g2_name, output_dir, result_filename)
    else:
        methods.with_PA_peaks(chromosomes, input1_dir, input2_dir, g1_name, g2_name, output_dir, result_filename)

# Print the total time taken and write the results to an Excel file
print("Total time:", round((time.time() - startTime) / 60, 2), "minutes.")
read_file = pd.read_csv(os.path.join(output_dir, result_filename + ".csv"), delimiter="\t")
read_file.to_excel(os.path.join(output_dir, result_filename + ".xlsx"), index=None, header=True)
