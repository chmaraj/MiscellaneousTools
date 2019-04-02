import pandas as pd
import numpy as np

"""Use this to convert files between tsv format and phylip"""


class Converter(object):

    def __init__(self, args):
        self.args = args
        self.in_file = args.input
        self.out_file = args.output

        if args.format == 'tsv':
            self.run_phylip_to_tsv()
        elif args.format == 'phylip':
            self.run_tsv_to_phylip()

    def run_phylip_to_tsv(self):
        with open(self.in_file, "r") as file:
            lines = file.read().strip().split("\n")  # Remove whitespace, and convert to list of lines
            sample_names = set()
            data = []
            for line in lines:
                if line == lines[0]:  # We don't need the first line, only contains number of samples
                    continue
                entries = line.split()
                name = entries[0]  # Name is the first field
                sample_names.add(name)
                data_row = entries[1:]  # Data contained in the remaining fields
                data.append(data_row)
            np_data = np.array(data)  # Convert data to numpy array so we can make a dataframe out of it
            sorted_names = sorted(sample_names)
            df = pd.DataFrame(np_data, index=sorted_names, columns=sorted_names)  # Make the dataframe
            df.to_csv(path_or_buf=self.out_file, sep="\t")  # Convert the dataframe to the desired tsv file

    def run_tsv_to_phylip(self):
        """Not implemented yet"""
        with open(self.in_file, "r") as file:
            lines = file.read().split("\n")


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Used to convert from .tsv to .phylip or vice-versa")
    parser.add_argument("-i", "--input", metavar="Input file", required=True,
                        help="Either a .tsv file or a .phylip file to convert")
    parser.add_argument("-o", "--output", metavar="Output file name", required=True,
                        help="Full desired path of final file")
    parser.add_argument("-f", "--format", metavar="tsv", required=True,
                        help="Desired output format. Currently only phylip_to_tsv is implemented,"
                             "so please use tsv.")

    arguments = parser.parse_args()

Converter(arguments)

