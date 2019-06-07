#!/usr/bin/env python
# coding: utf-8


from subprocess import Popen, PIPE
import os
import pathlib
import shutil
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genome", help="Genome fasta file", required=True)
parser.add_argument("-m", "--mites", help="MITEs fasta file", required=True)
parser.add_argument("-o", "--organism", help="Organism name", required=True)
args = parser.parse_args()

path_genome = args.genome
path_gff = args.mites
path_out_dir = 'data/' + args.organism + '/by_family/'

blastn -task blastn -query data/rice/2892.fa -db /home/juan/Desktop/juan/bio/data/genomes/rice/db -outfmt  "6 qseqid sseqid qstart qend sstart send mismatch gaps pident evalue length qlen slen qcovs" > data/rice/2892.genome.csv