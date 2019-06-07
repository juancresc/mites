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
parser.add_argument("-c", "--clusters", help="Clusters dir", required=True)
parser.add_argument("--csv", help="CSV clusters file", required=True)
parser.add_argument("--fasta", help="Fasta clusters file", required=True)
parser.add_argument("--data", help="data dir", required=True)
parser.add_argument("-g", "--genome", help="", required=True)
parser.add_argument("-i", "--pid", help="", default=0.97)
parser.add_argument("-m", "--min", help="Min elements in clusters", default=30)
parser.add_argument("-w", "--workers", help="Min elements in clusters", default=1)
args = parser.parse_args()

files = [f for f in os.listdir(args.families) if os.path.isfile(os.path.join(args.families, f))]
count = 0
total = len(files)
for file in files:
    file_path = args.families + file
    count += 1
    print("%i / %i" % (count,total))
    c_dir = args.clusters + file + "/"
    if os.path.exists(c_dir):
        continue
    if os.path.isdir(c_dir):
        shutil.rmtree(c_dir)
    pathlib.Path(c_dir).mkdir(parents=True, exist_ok=True)
    cmd_list = [
    './bin/vsearch-2.11.1/bin/vsearch',
    '--cluster_fast', args.output + k + ".flanking.fasta",
    '--threads',str(args.workers),
    '--strand','both',
    '--clusters', c_dir + "c_",
    '--iddef','1',
    '--id', str(args.pid)]
    p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
    out,err = p.communicate()




files = [f for f in os.listdir(args.families) if os.path.isfile(os.path.join(args.families, f))]
count = 0
total = len(files)
res = {}
valid_clusters = {}
for file in files:
    file_path = args.families + file
    c_dir = args.clusters + file + "/"
    clusters = [f for f in os.listdir(c_dir) if os.path.isfile(os.path.join(c_dir, f))]
    seqs = []
    for cluster in clusters:
        file_c = open(c_dir + cluster, "r")
        count_seqs = 0
        for line in file_c:
            if line[0:1] == ">":
                count_seqs += 1
        seqs.append(count_seqs)
        if count_seqs > args.min:
            valid_clusters[c_dir + cluster] = count_seqs
    if max(seqs) > args.min:
        res[file] = (len(clusters), max(seqs))
    count += 1
    print("%i / %i" % (count,total))

res_s = sorted(res.items(), key=lambda x:x[1][1])


with open(args.csv, 'w') as f:
    for v in res_s:
        f.write("%s,%s,%s\n"%(v[0],v[1][0],v[1][1]))
