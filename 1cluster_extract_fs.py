#!/usr/bin/env python
# coding: utf-8

import argparse
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
parser.add_argument("-f", "--families", help="Families dir", required=True)
parser.add_argument("--clusters", help="Clusters dir", required=True)
parser.add_argument("--fsdir", help="Flanking sequence dir", required=True)
parser.add_argument("--fslen", help="Flanking sequence length", default=50)
parser.add_argument("--with_fs", help="With fs", required=True)
parser.add_argument("--csv", help="CSV clusters file", required=True)
parser.add_argument("--fasta", help="Fasta clusters file", required=True)
parser.add_argument("--data", help="data dir", required=True)
parser.add_argument("-g", "--genome", help="", required=True)
parser.add_argument("-i", "--pid", help="", default=0.96)
parser.add_argument("-m", "--min", help="Min elements in clusters", default=30)
parser.add_argument("-w", "--workers", help="Min elements in clusters", default=1)
args = parser.parse_args()

if os.path.isdir(args.fsdir):
    shutil.rmtree(args.fsdir)
pathlib.Path(args.fsdir).mkdir(parents=True, exist_ok=True)

if os.path.isdir(args.with_fs):
    shutil.rmtree(args.with_fs)
pathlib.Path(args.with_fs).mkdir(parents=True, exist_ok=True)


df = pd.read_csv(args.csv, sep=',', header=None)
df.columns = ['family','file','cluster','counts']

fasta_seq = SeqIO.parse(args.genome, 'fasta')
buffer_mites_fs = {}
buffer_mites_ofs = {}
print("extracting sequences...")
for record in fasta_seq:
    for k,v in df.iterrows():
        buffer_with_fs = []
        buffer_fs = []
        fasta_seq = SeqIO.parse(v.file, 'fasta')
        for record in fasta_seq:
            seqname_data = record.id.split('_')
            chromosome = seqname_data[-3]
            start = int(seqname_data[-2])
            end = int(seqname_data[-1])
            if chromosome != record.id:
                continue
            clean_seq = ''.join(str(record.seq).splitlines())

            #add flanking
            start_f = max(start - args.fslen,0)
            end_f = end + args.fslen
            new_seq = clean_seq[start_f:end_f]
            id = record.id + "_fs_" + str(start) + '_' + str(end)
            seq = SeqRecord(Seq(new_seq), id=id, description="_")
            buffer_with_fs.append(seq)

            #only flanking 1
            start_o_f_1 = max(start - int(args.fslen),0)
            end_o_f_1 = start
            new_seq = clean_seq[start_o_f_1:end_o_f_1]
            id = record.id + "_flank_1_" + str(start_o_f_1) + '_' + str(end_o_f_1)
            seq = SeqRecord(Seq(new_seq), id=id, description="_")
            buffer_fs.append(seq)
            
            #only flanking 2
            start_o_f_2 = end 
            end_o_f_2 = end + int(args.fslen)
            new_seq = clean_seq[start_o_f_2:end_o_f_2]
            id = record.id + "_flank_2_" + str(start_o_f_2) + '_' + str(end_o_f_2)
            seq = SeqRecord(Seq(new_seq), id=id, description="_")
            buffer_fs.append(seq)

        SeqIO.write(buffer_fs, args.fsdir + v.family + '/' + v.cluster + ".fasta", "fasta")
        print('written',args.fsdir + v.file + '_' + cluster)

        SeqIO.write(buffer_with_fs, args.with_fs + v.family + '/' + v.cluster + ".fasta", "fasta")
        print('written',args.with_fs + cluster)