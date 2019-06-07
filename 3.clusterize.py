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
path_out_cluster = 'data/' + args.organism + '/by_cluster/'
pid = 0.97
min_count = 10
workers = 5
path_out_dir_clus = 'data/' + args.organism + '/clusters_similar/'
pid = 0.97
path_out_fasta = 'data/' + args.organism + '/similar.fasta'
path_similar = 'data/' + args.organism + '/similar_mites_count.csv'

files = [f for f in os.listdir(path_out_dir) if os.path.isfile(os.path.join(path_out_dir, f))]
count = 0
total = len(files)
for file in files:
    file_path = path_out_dir + file
    count += 1
    print("%i / %i" % (count,total))
    c_dir = path_out_cluster + file + "/"
    if os.path.exists(c_dir):
        continue
    if os.path.isdir(c_dir):
        shutil.rmtree(c_dir)
    pathlib.Path(c_dir).mkdir(parents=True, exist_ok=True)
    cmd_list = [
    './bin/vsearch-2.11.1/bin/vsearch',
    '--cluster_fast', file_path,
    '--threads',str(workers),
    '--strand','both',
    '--clusters', c_dir + "c_",
    '--iddef','1',
    '--id', str(pid)]
    p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
    out,err = p.communicate()




files = [f for f in os.listdir(path_out_dir) if os.path.isfile(os.path.join(path_out_dir, f))]
count = 0
total = len(files)
res = {}
for file in files:
    file_path = path_out_dir + file
    c_dir = path_out_cluster + file + "/"
    clusters = [f for f in os.listdir(c_dir) if os.path.isfile(os.path.join(c_dir, f))]
    seqs = []
    for cluster in clusters:
        file_c = open(c_dir + cluster, "r")
        count_seqs = 0
        for line in file_c:
            if line[0:1] == ">":
                count_seqs += 1
        seqs.append(count_seqs)
    if max(seqs) > min_count:
        res[file] = (len(clusters), max(seqs))
        #print(file + ": " + str(len(clusters)) + " max: " + str(max(seqs)))
        #print("*" * 10)
    count += 1
    print("%i / %i" % (count,total))



res_s = sorted(res.items(), key=lambda x:x[1][0])


with open(path_similar, 'w') as f:
    for v in res_s:
        f.write("%s,%s,%s\n"%(v[0],v[1][0],v[1][1]))

df_similar = pd.read_csv(path_similar, sep=',', header=None)

df_similar.columns = ['MITE','clusters','counts']

df_similar = df_similar[df_similar.counts > 100]

buffer_mites = []
for k,v in df_similar.iterrows():
    mite = v.MITE
    count = v.counts
    file_path = path_out_dir + mite
    fasta_seq = SeqIO.parse(file_path, 'fasta')
    first_record = next(fasta_seq)
    buffer_mites.append(first_record)
    print(mite, count)

SeqIO.write(buffer_mites, path_out_fasta, "fasta")


cmd_list = [
'./bin/vsearch-2.11.1/bin/vsearch',
'--cluster_fast', path_out_fasta,
'--threads',str(workers),
'--strand','both',
'--clusters', path_out_dir_clus + "c_",
'--iddef','1',
'--id', '0.6']
p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
out,err = p.communicate()
