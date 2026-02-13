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

# create dictionary of degenerate bases
deg_bases = {
    "T": "T|U",
    "U": "T|U",
    "R": "A|G",
    "Y": "C|T",
    "N": "A|C|T|G|U"
}


# CLASSES -------------------------------------------------------------------------------------------

class ListofGenes:
    def __init__(self):
        self.genes = []

    def set_up_gene(self, header, sequence):
        gene_obj = Gene(header, sequence)
        self.genes.append(gene_obj)
    
    def draw_gene_base(self, x_margin, gene_height):

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
            ctx.set_line_width(4)
            # round ends of line
            ctx.set_line_cap(cairo.LINE_CAP_ROUND)

            # calculate y-location of line
            y_loc = (ind * (gene_height * 3)) + (gene_height * (1.5))
            setattr(gene, 'y_loc', y_loc)

            # define location of line
            ctx.move_to(x_margin, y_loc) # start
            ctx.line_to(x_margin + gene.length, y_loc) # end

            # draw line
            ctx.stroke()

        # draw exons
        for gene in self.genes:

            ctx.rectangle(x_margin + gene.exon_start, gene.y_loc - (gene_height/2), gene.exon_len, gene_height)

            ctx.fill()
        

        return surface


        
class Gene:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        self.length = len(sequence)
        self.get_exon_loc()


    def get_exon_loc(self):
        exon_match = re.search("[A-Z]+", sequence)
        self.exon_start =  exon_match.start() + 1
        self.exon_len = exon_match.end() - exon_match.start()
        # captialize entire gene sequence
        #self.sequence = self.sequence.upper()



colors = [(1, 0.984, 0), (1, 0.604, 0.231), (1, 0.604, 0.98), (0.078 , 0.439, 0.42)]

class ListofMotifs:
    def __init__(self):
        self.motifs = []

    def set_up_motif(self, motif_seq):
        motif_obj = Motif(motif_seq)
        self.motifs.append(motif_obj)

    def find_motifs(self, gene_objs: list[Gene]):
        for motif in self.motifs:
            motif.find_motifs(gene_objs)

    def draw_motifs(self, surface, x_margin, gene_height):
        ctx = cairo.Context(surface)
        for ind, motif in enumerate(self.motifs):
            ctx.set_source_rgb(*colors[ind])
            print(*colors[ind])
            for gene_y_loc in motif.motif_instances:
                for motif_instance in motif.motif_instances[gene_y_loc]:
                    
                    ctx.rectangle(x_margin + motif_instance + 1, gene_y_loc - (gene_height/2), motif.length, gene_height)
                    print(ctx.get_source().get_rgba())
                    ctx.fill()
                    #ctx.new_path()
        return surface




class Motif:
    def __init__(self, motif_seq):
        self.orig_motif = motif_seq
        self.length = len(motif_seq)
        self.get_motif_poss()
        # key is y-loc, value is list of start positions (1-based)
        self.motif_instances = {}
    
    def get_motif_poss(self):
        motif_upper = self.orig_motif.upper()
        self.motif_regex = "".join([("(" + deg_bases[char]+ ")") if char in deg_bases else char for char in motif_upper])

    def find_motifs(self, gene_objs: list[Gene]):
        for obj in gene_objs:
            motif_instances = re.finditer(self.motif_regex, obj.sequence, flags = re.IGNORECASE)
            self.motif_instances[obj.y_loc] = [instance.start() for instance in motif_instances]
            






# MAIN CODE -------------------------------------------------------------------------------------------

# initialize ListofGenes object
all_genes = ListofGenes()

# parse input FASTA file & collect each gene's information into a Gene object
with open(args.fasta, 'r') as fasta_file:
    # initialize a sequence variable to hold the entire sequence of an entry
    sequence = ""
    # loop through lines of FASTA file
    for ind, line in enumerate(fasta_file):
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

# initialize ListofMotifs object
all_motifs = ListofMotifs()

# parse input motifs file & initialize a Motif object for each motif
with open(args.motifs) as motifs_file:
    for line in motifs_file:
        line = line.strip("\n")
        all_motifs.set_up_motif(line)



surface = all_genes.draw_gene_base(20, 30)
all_motifs.find_motifs(all_genes.genes)
surface = all_motifs.draw_motifs(surface, 20, 30)

surface.write_to_png("motif_mark_output.png")