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
parser.add_argument("--onlyflanking", help="", required=True)
parser.add_argument("--flankingcluster", help="", required=True)
parser.add_argument("-w", "--workers", help="Min elements in clusters", default=1)
args = parser.parse_args()

clusters = [f for f in os.listdir(args.onlyflanking) if os.path.isfile(os.path.join(args.onlyflanking, f))]


for cluster in clusters:
    c_dir = args.flankingcluster + cluster + '/'
    if os.path.isdir(c_dir):
        shutil.rmtree(c_dir)
    pathlib.Path(c_dir).mkdir(parents=True, exist_ok=True)
    cmd_list = [
    './bin/vsearch-2.13.4/bin/vsearch',
    '--cluster_fast', args.onlyflanking + cluster,
    '--threads',str(args.workers),
    '--strand','both',
    '--clusters', c_dir + "c_",
    '--iddef','1',
    '--id', '0.3']
    p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
    out,err = p.communicate()
    file_c = open(args.onlyflanking + cluster, "r")
    count_seqs = 0
    for line in file_c:
        if line[0:1] == ">":
            count_seqs += 1
    files = [f for f in os.listdir(c_dir) if os.path.isfile(os.path.join(c_dir, f))]
    value = len(files) * 100 / (count_seqs)
    if value > 90:
        print(cluster, count_seqs,len(files), value)