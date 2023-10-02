import numpy as np
import argparse
import pandas as pd
import os
from scipy.signal import find_peaks

def remove_pos_bias(pos, direct, d):
    """
    remove position bias in pair file
    """
    shift = np.zeros(len(pos))
    shift[np.argwhere(direct=='-')] -= d
    shift[np.argwhere(direct=='+')] += d
    return pos + shift

parser = argparse.ArgumentParser(
    description=""" """)
parser.add_argument(
    'pair_file', help='contact pair file in tsv format, the first six columns are <contact id> <chromosome 1> <position 1> <chromosome 2> <position 2> <strand 1> <strand 2>', type=str)
parser.add_argument(
    '--path', help='output directory, default = ./', default='./', type=str)

args = parser.parse_args()
df = pd.read_csv(
    args.pair_file, header=None, sep='\t', comment='#')
chrom = str(df.iloc[0,1])
df = df[(df.iloc[:,1]==chrom) & (df.iloc[:,3]==chrom)] # contact within the same chromosome
df = df[abs(df.iloc[:,4] - df.iloc[:,2])>=100]
all_pos = np.array(df.iloc[:, [2, 4]]).flatten()
all_direct = np.array(df.iloc[:, [5, 6]].astype('str')).flatten()
all_direct2 = all_direct[::2]+all_direct[1::2]
d = df.iloc[:, 4] - df.iloc[:, 2]

s1, _ = np.histogram(d[all_direct2 == '+-'], np.linspace(0, 1000, 100))
s2, _ = np.histogram(d[all_direct2 == '++'], np.linspace(0, 1000, 100))
s3, _ = np.histogram(d[all_direct2 == '-+'], np.linspace(0, 1000, 100))
s4, _ = np.histogram(d[all_direct2 == '--'], np.linspace(0, 1000, 100))
all_pos_unbiased = remove_pos_bias(all_pos, all_direct, abs(find_peaks(
    s4, prominence=1000)[0][0] - find_peaks(s3, prominence=1000)[0][0])*5)
# all_pos_unbiased = remove_pos_bias(all_pos, all_direct, 70)

starts = all_pos_unbiased-80
ends = all_pos_unbiased+80
df = pd.DataFrame(
    {'chr': chrom.replace('chr', ''), 'starts': starts.astype(int), 'ends': ends.astype(int)})
df = df[np.logical_and(df.starts>0, df.ends>0)]

df.to_csv(os.path.join(args.path, 'reads.bed'), header=None, index=None, sep='\t')
np.save(os.path.join(args.path, 'all_pos_unbiased'), all_pos_unbiased)
