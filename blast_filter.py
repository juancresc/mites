import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="", required=True)
parser.add_argument("-o", "--output", help="", required=True)

parser.add_argument("--min_len", help="",default=50)
parser.add_argument("--max_len", help="",default=False)
parser.add_argument("--min_distance", help="",default=5)
parser.add_argument("--min_q", help="",default=0.8)
parser.add_argument("--max_q", help="",default=1.2)
parser.add_argument("--min_pident", help="",default=85)
parser.add_argument("--min_qcov", help="",default=85)
args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t', header=None)
df.columns = ['qseqid','sseqid','qstart','qend','sstart','send','mismatch','gaps','pident','evalue','length','qlen','slen','qcovs']
print('initial:',len(df.index))
initial = len(df.index)

#filter by length
if args.min_len:
    df = df[df.qlen > args.min_len]
    print('Min len: ' + str(len(df.index)))
    min_length = str(len(df.index))

if args.max_len:
    df = df[df.qlen < args.max_len]
    print('Max len: ' + str(len(df.index)))    
    max_length = str(len(df.index))


#filter by query / subject length treshold
if args.min_q:
    df = df[((df.length / df.qlen) >= args.min_q)]
    print('min treshold:',len(df.index))
    min_treshold = str(len(df.index))

if args.max_q:
    df = df[((df.length / df.qlen) <= args.max_q)]
    print('max treshold:',len(df.index))
    max_treshold = str(len(df.index))

#filter by pident
if args.min_pident:
    df = df[(df.pident >= args.min_pident)]
    print('Min_pident: ' + str(len(df.index)))
    min_pident = str(len(df.index))


#filter by qcov
if args.min_qcov:
    df = df[(df.qcovs >= args.min_qcov)]
    print('Min qcov: ' + str(len(df.index)))
    min_qcov = str(len(df.index))

if args.min_distance:
    print('Filtering overlapped...')
    #order sstart and send
    df['new_sstart'] = df[['sstart','send']].min(axis=1)
    df['new_ssend'] = df[['sstart','send']].max(axis=1)
    df['sstart'] = df['new_sstart']
    df['send'] = df['new_ssend']
    df = df.drop('new_sstart',axis=1).drop('new_ssend',axis=1)
    df = df.sort_values(by=['sseqid','sstart', 'send'])
    df = df.reset_index(drop=True)
    # sep by chr
    dfs = {}
    for seq in df.sseqid.unique():
        dfs[seq] = df[df.sseqid == seq]

    # filter overlapped 
    rows = []
    discard = []
    total = len(df.index)
    count = 0
    curr = 0
    for index, row in df.iterrows():
        count += 1
        curr_new = int(count * 100 * 1.0 / (total * 1.0))
        if curr_new != curr:
            curr = curr_new
            if curr_new % 5 == 0:
                print(curr_new)
        if index in discard:
            continue
        for k2, v2 in df.loc[index:,].iterrows():
            if abs(v2.sstart - row.sstart) > args.min_distance:
                break
            if abs(v2.sstart - row.sstart) <= args.min_distance and abs(v2.send - row.send) <= args.min_distance:
                discard.append(k2)
        rows.append(row)


    df = pd.DataFrame(rows)
    print('Non overlapped: ' + str(len(df.index)))
    non_overlapped = str(len(df.index))


df.to_csv(args.output, index=None, sep='\t')

print('Initial: ' + str(initial))
if args.min_len:
    print('Min len: ' + str(min_length))
if args.max_len:
    print('Max len: ' + str(max_length))
if args.min_q:
    print('Min treshold: ' + str(min_treshold))
if args.max_q:
    print('Max treshold: ' + str(max_treshold))
if args.min_pident:
    print('Min pident: ' + str(min_pident))
if args.min_qcov:
    print('Min qcov: ' + str(min_qcov))
if args.min_distance:
    print('Non overlapped: ' + str(non_overlapped))
print('Saved to: ' + args.output)
