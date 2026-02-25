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



# DEFINING CLASSES -------------------------------------------------------------------------------------------

class ListofGenes:
    def __init__(self):
        self.genes: list[Gene] = []

    def set_up_gene(self, header, sequence):
        # create a Gene object
        gene_obj = Gene(header, sequence)
        # append object to genes attribute (so that each gene can be kept track of)
        self.genes.append(gene_obj)
    
    def draw_gene_base(self):
        # count number of genes
        num_genes = len(self.genes)
        # get length of longest gene
        longest_gene_len = max([gene.length for gene in self.genes])
        # width of surface = length of longest gene plus margin on each side
        surface_width = (longest_gene_len * 4) + (2 * x_margin)
        # height of surface = space + vertical padding for each gene + extra space at bottom + top for title & key
        surface_height = num_genes * (gene_height * 3) + 490

        # set up surface
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, surface_width, surface_height)
        ctx = cairo.Context(surface)
        # make background white
        ctx.set_source_rgb(1, 1, 1)
        ctx.rectangle(0, 0, surface_width, surface_height)
        ctx.fill()
        # set color to dark blue for drawing gene
        ctx.set_source_rgb(0.121, 0.165, 0.212)

        # draw introns
        for ind, gene in enumerate(self.genes):
            # set thickness of line
            ctx.set_line_width(6)
            # round ends of line
            ctx.set_line_cap(cairo.LINE_CAP_ROUND)
            # calculate y-location of line
            y_loc = (ind * (gene_height * 3)) + (gene_height * (1.5)) + 285
            setattr(gene, 'y_loc', y_loc)
            # define location of line
            ctx.move_to(x_margin, y_loc) # start
            ctx.line_to(x_margin + (gene.length * 4), y_loc) # end
            # draw line
            ctx.stroke()

        # draw exons
        for gene in self.genes:
            # define location + dimensions of exon
            ctx.rectangle(x_margin + (gene.exon_start * 4), gene.y_loc - (gene_height/2), gene.exon_len * 4, gene_height)
            # draw rectangle
            ctx.fill()

        return surface
    
    def draw_headers(self, surface):
        ctx = cairo.Context(surface)
        for gene in self.genes:
            # set color to black for header text
            ctx.set_source_rgb(0, 0, 0)
            # define font aesthetics for gene name
            ctx.select_font_face("monospace", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
            ctx.set_font_size(70)
            # move to start position of gene heading
            ctx.move_to(x_margin, gene.y_loc - gene_height)
            # specify gene name to write out (first word in header)
            ctx.show_text(gene.header.strip(">").split(" ")[0])
            # define font aesthetics for subheading
            ctx.select_font_face("monospace", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            ctx.set_font_size(40)
            # move to start position of subheading
            ctx.move_to(x_margin, gene.y_loc - gene_height + 40)
            # specify subheading to write out (everything after gene name)
            ctx.show_text(" ".join(gene.header.split(" ")[1:]))

        return surface


class Gene:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        # calculate length of gene
        self.length = len(sequence)
        self.get_exon_loc()

    def get_exon_loc(self):
        # find exon using regex to match uppercase sequence
        exon_match = re.search("[A-Z]+", sequence)
        # extract start position of exon (1-based)
        self.exon_start =  exon_match.start() + 1
        # calculate exon length using start and end indices
        self.exon_len = exon_match.end() - exon_match.start()


class ListofMotifs:
    def __init__(self):
        self.motifs = []

    def set_up_motif(self, motif_seq):
        # create a Motif object
        motif_obj = Motif(motif_seq)
        # append object to motifs attribute (so that each gene can be kept track of)
        self.motifs.append(motif_obj)

    def find_motifs(self, gene_objs: list[Gene]):
        for motif in self.motifs:
            # find start locations of motif occurrences in gene
            motif.find_motifs(gene_objs)
    
    def draw_motifs(self, surface):
        ctx = cairo.Context(surface)
        # loop through Motif objects
        for ind, motif in enumerate(self.motifs):
            # assign unique color to motif
            ctx.set_source_rgb(*colors[ind])
            # loop through genes
            for gene_y_loc in motif.motif_instances:
                # loop through motif start locations
                for motif_instance in motif.motif_instances[gene_y_loc]:
                    # calculate y location of motif
                    motif_y_loc = gene_y_loc - (gene_height/2) + ((ind) * gene_height/num_motifs)
                    # draw motif as a rectangle
                    ctx.rectangle(x_margin + (motif_instance + 1) * 4, motif_y_loc, motif.length * 4, gene_height/len(self.motifs))
                    ctx.fill()
                    # draw line that marks end of motif
                    ctx.set_line_width(2)
                    ctx.move_to(x_margin + (motif_instance + 1) * 4 + (motif.length * 4) - 1, motif_y_loc - 10)
                    ctx.line_to(x_margin + (motif_instance + 1) * 4 + (motif.length * 4) -1 , motif_y_loc + 10)
                    ctx.stroke()

        return surface
    
    def draw_key(self, surface, motif_height):
        ctx = cairo.Context(surface)
        # create a variable to keep track of horizontal space used up by key so far
        current_x_loc = x_margin
        # y location of key
        y_loc = surface.get_height() - 120
        # loop through motifs
        for ind, motif in enumerate(self.motifs):
            # draw motif box (to scale)
            ctx.set_source_rgb(*colors[ind])
            ctx.rectangle(current_x_loc, y_loc, motif.length * 4, motif_height)
            ctx.fill()
            # add line signifying end of motif (to match how motifs are drawn within figure)
            ctx.set_line_width(2)
            ctx.move_to(current_x_loc + (motif.length * 4) - 1, y_loc)
            ctx.line_to(current_x_loc + (motif.length * 4) - 1, y_loc - 10)
            ctx.stroke()
            current_x_loc += ((motif.length * 4) + 13)
            # write out motif sequence/name
            ctx.set_source_rgb(0, 0, 0)
            ctx.select_font_face("monospace", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            ctx.set_font_size(40)
            ctx.move_to(current_x_loc, y_loc + 35)
            ctx.show_text(motif.orig_motif)
            # calculate space used up by motif name
            _, _, text_width, _, _, _ = ctx.text_extents(motif.orig_motif)
            current_x_loc += (text_width + 110)
        # add bar separating key from genes
        ctx.set_source_rgb(0.121, 0.165, 0.212)
        ctx.rectangle(0, y_loc - 190, surface.get_width(), 125)
        ctx.fill()
        # add header for key
        ctx.set_source_rgb(1, 1, 1)
        ctx.select_font_face("monospace", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(70)
        ctx.move_to(x_margin, y_loc - 109)
        ctx.show_text("Motif Key")

        return surface
    

class Motif:
    def __init__(self, motif_seq):
        self.orig_motif = motif_seq
        self.length = len(motif_seq)
        self.get_motif_poss()
        # key is y-loc, value is list of start positions (1-based)
        self.motif_instances = {}
    
    def get_motif_poss(self):
        # make motif upper case
        motif_upper = self.orig_motif.upper()
        # create regex pattern for motif, containing all possible combinations of matching letters
        self.motif_regex = "(?=" + "".join([("(" + deg_bases[char]+ ")") if char in deg_bases else char for char in motif_upper]) + ")"

    def find_motifs(self, gene_objs: list[Gene]):
        # loop through each gene
        for obj in gene_objs:
            # find all (including overlapping) matches of motif
            motif_instances = re.finditer(self.motif_regex, obj.sequence, flags = re.IGNORECASE)
            # get start locations of motifs and store in motif_instances dictionary, where key is y-location of gene
            self.motif_instances[obj.y_loc] = [instance.start() for instance in motif_instances]
            


# SET UP / INITIALIZATION -------------------------------------------------------------------------------------------

# create dictionary of degenerate bases
deg_bases = {
    "T": "T|U",
    "U": "T|U",
    "R": "A|G",
    "Y": "C|T",
    "N": "A|C|T|G|U"
}

# define a function to draw header at top of surface
def draw_header(surface, fasta_file):
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(0.121, 0.165, 0.212)
    ctx.rectangle(0, 0, surface.get_width(), 250)
    ctx.fill()
    fasta_title = fasta_file.split("/")[-1]
    ctx.set_source_rgb(1, 1, 1)
    ctx.select_font_face("monospace", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_font_size(105)
    ctx.move_to(x_margin, 165)
    ctx.show_text(f"Motifs on Genes From {fasta_title}")
    return surface

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

# parse input motifs file 
with open(args.motifs) as motifs_file:
    for line in motifs_file:
        line = line.strip("\n")
        # initialize a Motif object for each motif
        all_motifs.set_up_motif(line)

# define colors for motifs
colors = [(0.921, 0.788, 0.725), (0.573, 0.678, 0.643), (0.867, 0.557, 0.345), (0.996, 0.847, 0.651), (0.784, 0.765, 0.855)]

# define sizes for surface
x_margin = 140
gene_height = 200

# extract output file prefix from input fasta file name by removing file extension
out_file_prefix = args.fasta.split(".")[-2]


# MAIN CODE: CREATE DRAWING -------------------------------------------------------------------------------------------

# create surface + draw exons + introns
surface = all_genes.draw_gene_base()
# add header bar to top of surface
surface = draw_header(surface, args.fasta)
# find + store all instances of motifs in genes
all_motifs.find_motifs(all_genes.genes)
# count number of motifs
num_motifs = len(all_motifs.motifs)
# draw motifs
surface = all_motifs.draw_motifs(surface)
# add key for motifs at bottom of surface
surface = all_motifs.draw_key(surface, gene_height/num_motifs)
# add headers for each gene to surface
surface = all_genes.draw_headers(surface)


# CREATING OUTPUT FILES -------------------------------------------------------------------------------------------

# save drawing as output png file
surface.write_to_png(f"{out_file_prefix}.png")

# write motif counts per gene to a markdown file
with open(f"{out_file_prefix}_stats.md", 'w') as fh:
    # create header line
    fh.write("| |")
    for motif_seq in [motif_obj.orig_motif for motif_obj in all_motifs.motifs]:
        fh.write(f"{motif_seq}|")
    fh.write("\n|")
    # create delimiter row
    for motif in range(num_motifs + 1):
        fh.write("-|")
    fh.write("\n|")
    # loop through genes
    for ind, gene in enumerate(all_genes.genes):
        # extract gene name from header
        fh.write(f"{gene.header.strip(">").split(" ")[0]}|")
        # loop through motifs
        for motif in all_motifs.motifs:
            # count occurence of motif within gene
            fh.write(f"{len(motif.motif_instances[gene.y_loc])}|")
        fh.write("\n")
        # don't add new pipe character if last line of table has been written
        if ind < (len(all_genes.genes) - 1):
            fh.write("|")
