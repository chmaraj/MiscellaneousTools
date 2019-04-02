import subprocess
import os
import pathlib
from multiprocessing import Pool
from operator import itemgetter

base_dir = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CloselyRelatedSamples/Plasmid_High_Throughput"
# assemblies_dir = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CampSalm2/CampSalm2/Assemblies/" \
#                 "Salmonella_Similarities/Copy_of_Assemblies"
deltas_dir = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CloselyRelatedSamples/Plasmid_High_Throughput/Deltas"
alignments_dir = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CloselyRelatedSamples/Plasmid_High_Throughput/Alignments"
dna_diff_dir = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CloselyRelatedSamples/Plasmid_High_Throughput/DNADiff"
# dna_diff_dir = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CampSalm2/CampSalm2/Assemblies/Salmonella_Similarities/" \
#                "Genome/DNADiff"
plasmid_dir = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CloselyRelatedSamples/Plasmid_High_Throughput/Plasmid_Sequences"
related_dir = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CloselyRelatedSamples/Plasmid_High_Throughput/" \
              "ReferenceSamples/fna"

# sample_list = []
# for entry in os.listdir(assemblies_dir):
#        name = entry.split(".")[0]
#        sample_list.append(name)

# plasmid_list = []
# related_list = []
# for entry in os.listdir(plasmid_dir):
#     name = "".join(entry.split(".")[:-1])
#     plasmid_list.append(name)

# for entry in os.listdir(related_dir):
#     name = "".join(entry.split(".")[:-1])
#     related_list.append(name)


def nucmer_run(entry):
    ref = "{}/{}.fasta".format(assemblies_dir, entry)
    current_dir = "{}/{}".format(deltas_dir, entry)
    pathlib.Path(current_dir).mkdir(parents=True, exist_ok=True)
    non_ref_list = [item for item in sample_list if item != entry]
    for item in non_ref_list:
        query = "{}/{}.fasta".format(assemblies_dir, item)
        delta = "{}/{}.delta".format(current_dir, item)
        command = ["nucmer", "--mum", "--threads=12", "--delta={}".format(delta), ref, query]
        subprocess.run(command)


def mummerplot_run(entry):
    current_dir = "{}/Alignment_to_{}".format(alignments_dir, entry)
    pathlib.Path(current_dir).mkdir(parents=True, exist_ok=True)
    non_ref_list = [item for item in sample_list if item != entry]
    for item in non_ref_list:
        p = "{}/{}".format(current_dir, item)
        input_item = "{}/{}/{}.delta".format(deltas_dir, entry, item)
        command = ["mummerplot", "-p", p, "-title", "Alignment_Image", "--layout", "--large", "--color",
                   "--png", input_item]
        subprocess.run(command)


def dna_diff_run(entry):
    current_dir = "{}/DNADiff_from_{}".format(dna_diff_dir, entry)
    pathlib.Path(current_dir).mkdir(parents=True, exist_ok=True)
    non_ref_list = [item for item in sample_list if item != entry]
    for item in non_ref_list:
        p = "{}/{}".format(current_dir, item)
        d = "{}/{}/{}.delta".format(deltas_dir, entry, item)
        command = ["dnadiff", "-p", p, "-d", d]
        subprocess.run(command)


def extract_plasmid(entry):
    entry_file = "{}/{}.fasta".format(assemblies_dir, entry)
    plasmid_file = "{}/{}_plasmid.fasta".format(plasmid_dir, entry)
    plasmid = ""
    with open(entry_file, "r") as item:
        lines = item.read()
        split_lines = lines.split(">")
        for sequence in split_lines[1:]:
            header = sequence.split("\n")[0]
            print(header)
            length_entry = header.split()[1]
            length = length_entry.split("=")[1]
            if 86000 < int(length) < 87000:
                plasmid = sequence
    with open(plasmid_file, "w") as output:
        output.write(plasmid)


def nucmer_plasmid(entry):
    ref = "{}/{}.fasta".format(plasmid_dir, entry)
    current_dir = "{}/{}".format(deltas_dir, entry)
    pathlib.Path(current_dir).mkdir(parents=True, exist_ok=True)
    for item in related_list:
        print("Running nucmer on {}:{}".format(entry, item))
        query = "{}/{}.fna".format(related_dir, item)
        delta = "{}/{}.delta".format(current_dir, item)
        command = ["nucmer", "--mum", "--threads=12", "--delta={}".format(delta), ref, query]
        subprocess.run(command)


def mummerplot_plasmid(entry):
    current_dir = "{}/Alignment_to_{}".format(alignments_dir, entry)
    pathlib.Path(current_dir).mkdir(parents=True, exist_ok=True)
    for item in related_list:
        p = "{}/{}".format(current_dir, item)
        input_item = "{}/{}/{}.delta".format(deltas_dir, entry, item)
        command = ["mummerplot", "-p", p, "-title", "Alignment_Image", "--layout", "--large", "--color",
                   "--png", input_item]
        subprocess.run(command)

def dna_diff_plasmid(entry):
    current_dir = "{}/DNADiff_from_{}".format(dna_diff_dir, entry)
    pathlib.Path(current_dir).mkdir(parents=True, exist_ok=True)
    for item in related_list:
        p = "{}/{}".format(current_dir, item)
        d = "{}/{}/{}.delta".format(deltas_dir, entry, item)
        command = ["dnadiff", "-p", p, "-d", d]
        subprocess.run(command)


# pool = Pool(1)
# pool.map_async(nucmer_plasmid, plasmid_list)
# pool.close()
# pool.join()

# pool2 = Pool(12)
# pool2.map_async(mummerplot_plasmid, plasmid_list)
# pool2.close()
# pool2.join()

# pool3 = Pool(12)
# pool3.map_async(dna_diff_plasmid, plasmid_list)
# pool3.close()
# pool3.join()

def compile_reports(searched_directory):
    with open("{}/{}/compiled_reports_60.txt".format(dna_diff_dir, searched_directory), "w") as log:
        for entry in os.listdir("{}/{}".format(dna_diff_dir, searched_directory)):
            if ".report" in entry:
                with open("{}/{}/{}".format(dna_diff_dir, searched_directory, entry)) as current_file:
                    lines = current_file.read()
                    split_lines = lines.split("\n")
                    for line in split_lines:
                        if "UnalignedBases" in line:
                            split_line = line.split()
                            reference_entry = split_line[1]
                            percentage = reference_entry.split("(")[1].replace("%)", "")
                            if float(percentage) < 40:
                                log.write("{}/{}: {}".format(searched_directory, entry, percentage))
                                log.write("\n")

def compile_compiled_reports(searched_directory):
    with open("{}/{}/compiled_reports.txt".format(dna_diff_dir, searched_directory), "w") as log:
        percentiles = []
        for entry in os.listdir("{}/{}".format(dna_diff_dir, searched_directory)):
            if "compiled_reports_" in entry:
                with open("{}/{}/{}".format(dna_diff_dir, searched_directory, entry)) as current_file:
                    lines = current_file.read()
                    split_lines = lines.split("\n")
                    percentile = entry.split(".")[0].split("_")[-1]
                    percentiles.append((percentile, len(split_lines)))
        percentiles.sort(key=itemgetter(0), reverse=True)
        for entry in percentiles:
            log.write("{}th percentile: {} entries".format(entry[0], entry[1]))
            log.write("\n")


pool = Pool(6)
pool.map_async(compile_compiled_reports, os.listdir(dna_diff_dir))
pool.close()
pool.join()



