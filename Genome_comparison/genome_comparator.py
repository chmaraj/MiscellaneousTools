import os
from multiprocessing import Pool
import pathlib
import subprocess
import re
import shutil

import pandas as pd
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import linkage

"""A compilation of scripts to use with mash to compare genomes"""

base_dir = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/SalmonellaMASH/BacSort_subset_results"
assemblies_folder = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/SalmonellaMASH/BacSort_subset_results/Salmonella_001_cluster"
out_directory = "{}/mash/sketches".format(base_dir)
sketches_folder = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/SalmonellaMASH/BacSort_subset_results/mash/sketches"
pathlib.Path(out_directory).mkdir(parents=True, exist_ok=True)


def sketch_it(entry_name):
    # produce a sketch file for the input file
    in_file = "{}/{}".format(assemblies_folder, entry_name)
    out_file = "{}/{}.msh".format(out_directory, entry_name)
    command = ["mash", "sketch", "-k", "21", "-s", "10000", "-p", "12", "-o", out_file, in_file]
    subprocess.run(command)


def paste_it(sketch_list):
    # paste all sketch files together. Must be all contained on a text file as a list.
    out_file = "{}/all_sketches.msh".format(base_dir)
    command = ["mash", "paste", out_file, "-l", sketch_list]
    subprocess.run(command)


def dist_it(entry_name):
    # produce distances from each sample to all references
    in_file = "{}/mash/sketches/{}".format(base_dir, entry_name)
    reference_file = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/SalmonellaMASH/BacSort_subset_results/all_sketches.msh"
    command = ["mash", "dist", "-p", "12", "-t", reference_file, in_file]
    p1 = subprocess.Popen(command, stdout = subprocess.PIPE)
    with open("{}/mash/distances/{}_dist.tsv".format(base_dir, entry_name), "wb") as output:
        output.write(p1.stdout.read())
        p1.wait()


def trim_tsv(tsv_file):
    # trim extra information from samples names in the tsv files
    old_file = '/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/SalmonellaMASH/BacSort_subset_results/Salmonella_001_cluster/'
    with open("{}/mash/distances/{}".format(base_dir, tsv_file), "r") as entry:
        entry_contents = entry.read()
        new_contents = re.sub(old_file, '', entry_contents)
        new_contents2 = re.sub(r'.fna.gz', '', new_contents)
        new_contents3 = re.sub(r'#query', '', new_contents2)
        new_contents4 = re.sub("{}/".format(assemblies_folder), '', new_contents3)
        with open("{}/mash/distances/{}.tmp".format(base_dir, tsv_file), "w") as output:
            output.write(new_contents4)
    shutil.move("{}/mash/distances/{}.tmp".format(base_dir, tsv_file), "{}/mash/distances/{}".format(base_dir, tsv_file))


def rename_tsvs(directory):
    # rename tsvs if necessary
    tmp_files = [x for x in os.listdir(directory) if ".tmp" in x]
    for entry in tmp_files:
        shutil.move("{}/mash/distances/{}".format(base_dir, entry), "{}/mash/distances/{}".format(base_dir, entry.split(".")[:-1]))


def remove_files(directory, pattern):
    # remove excess files
    for entry in os.listdir(directory):
        if pattern in entry:
            os.remove("{}/{}".format(directory, entry))            


def concatenate_distances(directory):
    # create a tsv file containing all distances, generating a distance matrix for all samples
    with open("{}/all_dist.tsv".format(base_dir), "a") as output:
        with open("/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/SalmonellaMASH/BacSort_subset_results/mash/distances/06D1274-20-15.fna.gz.msh_dist.tsv", "r") as entry:
            output.write(entry.readline())
        for entry in os.listdir(directory):
            with open("{}/{}".format(directory, entry), "r") as entry:
                output.write(entry.read().split("\n")[1])
                output.write("\n")


def dendrogram_constructor(distance_matrix):
    # generate the hc.nwk tree
    dataframe = pd.read_csv(distance_matrix, sep='\t', index_col=0)
    print("Made dataframe")
    df = pd.DataFrame(index=dataframe.index)
    for entry in dataframe.index:
        df[entry] = dataframe[entry]
    np_df = df.to_numpy()
    print("Converted to np")
    linkage_matrix = linkage(np_df, "single")
    print("Completed linkage")
    dendrogram(linkage_matrix, color_threshold = 1, show_leaf_counts = True)
    plt.title = "Test"
    plt.show()


# pool = Pool(1)
# pool.map_async(dist_it, os.listdir(sketches_folder))
# pool.close()
# pool.join()

# dendrogram_constructor("/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CampSalmMASH/all_dist.tsv")

# pool = Pool(12)
# output = pool.map_async(trim_tsv, os.listdir("{}/mash/distances".format(base_dir)))
# pool.close()
# pool.join()

concatenate_distances("/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/SalmonellaMASH/BacSort_subset_results/mash/distances")


# with open("/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CampSalmMASH/plasmid_sketch.list", "w") as output:
#    for entry in os.listdir("/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CampSalmMASH/mash/all_sketches/Plasmid"):
#        output.write("/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/CampSalmMASH/mash/all_sketches/Plasmid/{}".format(entry))
#        output.write("\n")

# sketch_list = "/media/bioinfo3/3tb_hdd/ChmaraJ/CampSalm/SalmonellaMASH/BacSort_subset_results/all_sketches.txt"
# paste_it(sketch_list)

        
