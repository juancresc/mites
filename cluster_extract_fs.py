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

#list families and cluster them
files = [f for f in os.listdir(args.families) if os.path.isfile(os.path.join(args.families, f))]
count = 0
total = len(files)
for file in files:
    file_path = args.families + file
    count += 1
    print("%i / %i" % (count,total))
    c_dir = args.clusters + file + "/"
    if os.path.isdir(c_dir):
        continue
        shutil.rmtree(c_dir)
    pathlib.Path(c_dir).mkdir(parents=True, exist_ok=True)
    cmd_list = [
    './bin/vsearch-2.13.4/bin/vsearch',
    '--cluster_fast', file_path,
    '--threads',str(args.workers),
    '--strand','both',
    '--clusters', c_dir + "c_",
    '--iddef','1',
    '--id', str(args.pid)]
    p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
    out,err = p.communicate()

#list cluster and filter by min number
files = [f for f in os.listdir(args.families) if os.path.isfile(os.path.join(args.families, f))]
count = 0
total = len(files)
res = {}
valid_clusters = {}
f = open(args.csv, 'w')
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
            #extract FS
            fasta_seq = SeqIO.parse(c_dir + cluster, 'fasta')
            buffer_fs = []
            buffer_with_fs = []
            for record in fasta_seq:
                seqname_data = record.id.split('_')
                chromosome = seqname_data[-3]
                start = int(seqname_data[-2])
                end = int(seqname_data[-1])

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

            
            SeqIO.write(buffer_fs, args.fsdir + file + '_' + cluster + ".fasta", "fasta")
            print('written',args.fsdir + file + '_' + cluster)

            SeqIO.write(buffer_with_fs, args.with_fs + file + '_' + cluster + ".fasta", "fasta")
            print('written',args.with_fs + cluster)
            #save results in a csv
            clus_file = c_dir + cluster
            f.write("%s,%s,%s\n"%(file,clus_file,count_seqs))
            
#            print(c_dir + cluster, count_seqs)
#    if len(seqs) > 0:
#        if max(seqs) > args.min:
#            res[file] = (len(clusters), max(seqs))
#        count += 1
#        print("%i / %i" % (count,total))

#res_s = sorted(res.items(), key=lambda x:x[1][1])


#with open(args.csv, 'w') as f:
#    for v in res_s:
#        f.write("%s,%s,%s\n"%(v[0],v[1][0],v[1][1]))
