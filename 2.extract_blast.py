#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

df = pd.read_csv(path_gff, sep='\t')

df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
print(len(df.index))


fasta_seq = SeqIO.parse(path_genome, 'fasta')
buffer_mites = {}

for record in fasta_seq:
    dff_extract = df[df.seqname == record.id]
    count = 0
    curr = 0
    total = len(dff_extract.index)
    print(record.id, len(dff_extract.index))
    clean_seq = ''.join(str(record.seq).splitlines())
    for k,v in dff_extract.iterrows():
        start = min(v.start,v.end)
        end = max(v.start,v.end)
        new_seq = clean_seq[start:end]
        att = v.attribute
        id = att + "_" + record.id + "_" + str(start) + '_' + str(end)
        seq = SeqRecord(Seq(new_seq), id=id, description="_")
        if not att in buffer_mites:
            buffer_mites[att] = []
        buffer_mites[att].append(seq)
        count += 1
        curr_new = int(count * 100 * 1.0 / (total * 1.0))
        if curr_new != curr:
            curr = curr_new
            if curr_new % 1 == 0:
                print(curr_new)

#writes fasta of all mites found in the gff
for k,v in buffer_mites.items():
    SeqIO.write(v, path_out_dir + k + ".fasta", "fasta")

