#!/usr/bin/env python3

# import libraries
import argparse
import cairo

# set global variables to hold inputs
def get_args():
    parser = argparse.ArgumentParser(description="A script that visualizes motifs on genes.")
    parser.add_argument("-f", "--fasta", help="input FASTA filename", required=True, type=str)
    parser.add_argument("-m", "--motifs", help="input motifs filename", required=True, type=str)
    return parser.parse_args()
args = get_args()
print(args)


# CLASSES -------------------------------------------------------------------------------------------

class ListofGenes:
    def __init__(self, list_of_genes: list[Gene]):
        self.genes = list_of_genes

    def set_up_gene(self, header, sequence):
        gene_obj = Gene(header, sequence)
        self.genes.append(gene_obj)

        
class Gene:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence





# MAIN CODE -------------------------------------------------------------------------------------------

# initialize ListofGenes object
all_genes = ListofGenes([])


# parse input FASTA file & collect each gene's information into a Gene object
with open(args.fasta, 'r') as fasta:
    # initialize a sequence variable to hold the entire sequence of an entry
    sequence = ""
    # loop through lines of FASTA file
    for ind, line in enumerate(fasta):
        line = line.strip("\n")
        # save first line as header
        if ind == 0:
            header = line
            continue
        # reaching header signifies that the entire sequence for the previous FASTA record has been collected
        if ">" in line:
            # set up Gene object for entry
            all_genes.set_up_gene(header, sequence)
            # re-initialize sequence variable to hold entire sequence of next entry
            sequence = ""
            # overwrite header to be next entry's header
            header = line
        else:
            # append the sequence to the existing sequence
            sequence += line
    # set up Gene object for last entry in FASTA file
    all_genes.set_up_gene(header, sequence)