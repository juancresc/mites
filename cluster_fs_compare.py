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
parser.add_argument("--csv", help="CSV clusters file", required=True)
parser.add_argument("-f", "--families", help="Dir", required=True)
parser.add_argument("-c", "--fsclusters", help="", required=True)
parser.add_argument("-w", "--workers", help="Min elements in clusters", default=1)
args = parser.parse_args()

df = pd.read_csv(args.csv, sep=',', header=None)
df.columns = ['MITE','clusters','counts']

buffer_mites = []
for k,v in df.iterrows():
    mite = v.MITE
    c_dir = args.fsclusters + mite
    if os.path.isdir(c_dir):
        shutil.rmtree(c_dir)
    pathlib.Path(c_dir).mkdir(parents=True, exist_ok=True)
    print(args.families + mite + ".onlyflanking.fasta")
    cmd_list = [
    './bin/vsearch-2.11.1/bin/vsearch',
    '--cluster_fast', args.families + mite + ".onlyflanking.fasta",
    '--threads',str(args.workers),
    '--strand','both',
    '--clusters', c_dir + "c_",
    '--iddef','1',
    '--id', '0.3']
    p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
    out,err = p.communicate()

    files = [f for f in os.listdir(c_dir) if os.path.isfile(os.path.join(c_dir, f))]
    value = (len(files) / 2) * 100 / v.counts
    print(mite, value)