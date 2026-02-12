#!/usr/bin/env python3

# import libraries
import argparse
import re
import cairo


# set global variables to hold inputs
def get_args():
    parser = argparse.ArgumentParser(description="A script that visualizes motifs on genes.")
    parser.add_argument("-f", "--fasta", help="input FASTA filename", required=True, type=str)
    parser.add_argument("-m", "--motifs", help="input motifs filename", required=True, type=str)
    return parser.parse_args()
args = get_args()


# CLASSES -------------------------------------------------------------------------------------------

class ListofGenes:
    def __init__(self, list_of_genes: list[Gene]):
        self.genes = list_of_genes

    def set_up_gene(self, header, sequence):
        gene_obj = Gene(header, sequence)
        self.genes.append(gene_obj)
    
    def draw_gene_bases(self, x_margin, gene_height):

        num_genes = len(self.genes)
        longest_gene_len = max([gene.length for gene in self.genes])

        surface_width = longest_gene_len + (2 * x_margin)
        # extra space at bottom for key
        surface_height = (num_genes + 1) * (gene_height * 3)

        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, surface_width, surface_height)
        ctx = cairo.Context(surface)

        # black
        ctx.set_source_rgb(0, 0, 0)

        # draw introns
        for ind, gene in enumerate(self.genes):
        
            # set thickness of line
            ctx.set_line_width(2)
            # round ends of line
            ctx.set_line_cap(cairo.LINE_CAP_ROUND)

            # calculate y-location of line
            y_loc = (ind * (gene_height * 3)) + (gene_height * (1.5))

            # define location of line
            ctx.move_to(x_margin, y_loc) # start
            ctx.line_to(x_margin + gene.length, y_loc) # end

            # draw line
            ctx.stroke()

        return surface


        
class Gene:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        self.length = len(sequence)
        self.get_exon_loc()


    def get_exon_loc(self):
        exon_match = re.search(r"[A-Z]+", sequence)
        self.exon =  [exon_match.start() + 1, exon_match.end()]
        self.exon_seq = sequence[exon_match.start(): exon_match.end()]





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

    #print([gene.exon for gene in all_genes.genes])
    #print([gene.exon_seq for gene in all_genes.genes])


    surface = all_genes.draw_gene_bases(20, 20)
    surface.write_to_png("motif_mark_output.png")