#!/usr/bin/env python3

# import libraries
import argparse
import cairo

# set global variables to hold inputs
def get_args():
    parser = argparse.ArgumentParser(description="A script that visualizes motifs on genes.")
    parser.add_argument("-f", "fasta", help="input FASTA filename", required=True, type=str)
    parser.add_argument("-m", "motifs", help="input motifs filename", required=True, type=str)
    return parser.parse_args()
args = get_args()