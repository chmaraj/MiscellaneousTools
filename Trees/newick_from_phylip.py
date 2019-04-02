import numpy as np
import pandas as pd
import os
from os.path import basename
from time import time

"""Use this tool to convert a phylip file to a newick tree."""


class NwkfromPhylip(object):

    def __init__(self, args):
        self.args = args
        self.input_file = args.input
        self.output = args.output

        self.run()

    def run(self):
        print("Parsing input file to pandas dataframe...", end="", flush=True)
        with open(self.input_file, "r") as file:
            lines = file.read().strip()  # Strip end lines if necessary
            useful_lines = lines.split("\n")[1:]  # First line is just the number of entries
            cluster_names = set()
            data = []
            for line in useful_lines:
                name = line.split("\t")[0]
                cluster_names.add(name)  # This gives us the column and index names we need
                data_row = []
                for data_point in line.split('\t')[1:]:
                    data_row.append(data_point)
                data.append(data_row)  # We now have the data in a nested list
            numpy_data = np.array(data)  # Convert nested list to numpy array so we can turn it into a dataframe
            sorted_names = sorted(cluster_names)  # Sort names to be in correct order
            df = pd.DataFrame(numpy_data, index=sorted_names, columns=sorted_names)  # Make the dataframe
            self.make_tree(df)


    def make_tree(self, df):

        labels = df.columns.tolist()
        labels[:] = ['\'{}\''.format(x) for x in labels]  # Ensure labels are correct

        # Using scipy -> Takes about 8s to run with ~9,000 Salmonella genomes (fasta)
        self.make_hc_dendrogram(df, labels)


    def make_hc_dendrogram(self, df, labels):
        from scipy.cluster.hierarchy import linkage, to_tree
        from scipy.spatial.distance import squareform

        # Convert square matrix to condensed distance matrix
        dists = squareform(df, checks=False)

        # hierarchical clustering
        linkage_matrix = linkage(dists, "ward")  # optimal_ordering=Ture is slow on large datasets

        # Create tree in Newick format
        print('Making dendrogram...', end="", flush=True)
        tree = to_tree(linkage_matrix, False)
        nw = self.getNewick(tree, "", tree.dist, labels)

        # save tree to file
        name = basename(self.input_file).split('.')[0]  # input file name without extension
        out_tree_file = self.output + '/' + name + '_hc.nwk'
        with open(out_tree_file, 'w') as out:
            out.write(nw)

    def getNewick(self, node, newick, parentdist, leaf_names):
        """
        https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
        :param node:
        :param newick:
        :param parentdist:
        :param leaf_names:
        :return:
        """
        if node.is_leaf():
            return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"

            newick = self.getNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = self.getNewick(node.get_right(), ",%s" % newick, node.dist, leaf_names)
            newick = "(%s" % newick

            return newick


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Convert a PHYLIP distance matrix into a Newick tree file')
    parser.add_argument('-i', '--input', metavar='input_file.phylip', required=True,
                        help='The input file to convert to a Newick tree')
    parser.add_argument('-o', '--output', metavar='Out folder', required=True,
                        help='Folder to hold the result file')

    arguments = parser.parse_args()

NwkfromPhylip(arguments)

