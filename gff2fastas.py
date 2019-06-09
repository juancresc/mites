import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="GFF file", required=True)
parser.add_argument("-g", "--genome", help="", required=True)
parser.add_argument("-o", "--output", help="Dir", required=True)
args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t', header=None)
df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

fasta_seq = SeqIO.parse(args.genome, 'fasta')
buffer_mites = {}
buffer_mites_fs = {}
buffer_mites_ofs = {}
print("extracting sequences...")
for record in fasta_seq:
    dff_extract = df[df.seqname == record.id]
    count = 0
    curr = 0
    total = len(dff_extract.index)
    print(record.id, len(dff_extract.index))
    clean_seq = ''.join(str(record.seq).splitlines())
    for k,v in dff_extract.iterrows():
        
        start = min(int(v.start),int(v.end))
        end = max(int(v.start),int(v.end))
        new_seq = clean_seq[start:end]
        att = v.attribute
        id = att + "_" + record.id + "_" + str(start) + '_' + str(end)
        seq = SeqRecord(Seq(new_seq), id=id, description="_")
        if not att in buffer_mites:
            buffer_mites[att] = []
        buffer_mites[att].append(seq)
        continue
        #add flanking
        start_f = max(start - args.flank_len,0)
        end_f = end + args.flank_len
        new_seq = clean_seq[start_f:end_f]
        att = v.attribute
        id = att + "_" + record.id + "_" + str(start) + '_' + str(end)
        seq = SeqRecord(Seq(new_seq), id=id, description="_")
        if not att in buffer_mites_fs:
            buffer_mites_fs[att] = []
        buffer_mites_fs[att].append(seq)
        
        #only flanking 1
        start_o_f_1 = max(start-args.flank_len,0)
        end_o_f_1 = start
        new_seq = clean_seq[start_o_f_1:end_o_f_1]
        att = v.attribute
        id = att + "_fs1_" + record.id + "_" + str(start) + '_' + str(end)
        seq = SeqRecord(Seq(new_seq), id=id, description="_")
        if not att in buffer_mites_ofs:
            buffer_mites_ofs[att] = []
        buffer_mites_ofs[att].append(seq)
        
        #only flanking 2
        start_o_f_2 = end 
        end_o_f_2 = end+args.flank_len
        new_seq = clean_seq[start_o_f_2:end_o_f_2]
        att = v.attribute
        id = att + "_fs2_" + record.id + "_" + str(start) + '_' + str(end)
        seq = SeqRecord(Seq(new_seq), id=id, description="_")
        buffer_mites_ofs[att].append(seq)
        
        count += 1
        curr_new = int(count * 100 * 1.0 / (total * 1.0))
        if curr_new != curr:
            curr = curr_new
            if curr_new % 10 == 0:
                print(curr_new)

for k,v in buffer_mites.items():
    SeqIO.write(v, args.output + k + ".fasta", "fasta")