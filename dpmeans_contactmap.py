import os
import numpy as np
from sys import argv
import argparse
from utils import dpmeans, encode_seq, fasta_iter

parser = argparse.ArgumentParser(
    description=""" """)
parser.add_argument('read_centers', help='read centers matrix in .npy format', type=str)
parser.add_argument(
    'path', help='output directory, default = ./', default='./', type=str)
parser.add_argument(
    'pwm', help='PWM of dinucleotide motifs, default = None', default=None, type=str)
parser.add_argument(
    'seq_file', help='reference genome fasta file, default = None', default=None, type=str)

args = parser.parse_args()
path = args.path
all_pos_unbiased = np.load(args.read_centers)
if args.pwm:
    pwm = np.load(args.pwm)
    seq_file = args.seqfile
part_pos = args.read_centers

bin_len = 1000
bins = np.arange(min(part_pos), max(part_pos)+bin_len, bin_len)
min_ends = 10
assignment = np.digitize(part_pos, bins=bins)-1
contact_assignment = (np.digitize(all_pos_unbiased, bins=bins)-1).reshape([-1,2])
c = all_pos_unbiased.reshape([-1,2])
n_peaks = 0
w = 250000 # cl element penalty
w2 = 10 #seq coef
peaks = []
pi = np.zeros(len(part_pos))
lbd = int(argv[5])
fiter = fasta_iter(seq_file)

for ff in fiter:
    header, seq = ff
    e = encode_seq(seq)
    e2 = encode_seq(seq[1:])
    scores = np.correlate(e, pwm[0], 'same') + np.correlate(1-e, pwm[1], 'same')
    scores2 = np.correlate(e2, pwm[0], 'same') + np.correlate(1-e2, pwm[1], 'same')
    scores = np.vstack([scores, np.r_[scores2, 0]]).T.flatten()
    scores = w2*(scores-scores.mean())/np.std(scores)

for i in range(len(bins)-1):
    bin_start = bins[i]
    if (i+1)%100 == 0:
        print(f'iteration: {i}, n_peaks: {n_peaks}, len(peaks): {len(np.concatenate(peaks))}, len(pi): {max(pi)+1}')
    idx = assignment==i
    ends = part_pos[idx]

    idx2 = (contact_assignment==i).all(1)
    cl_idx = c[idx2]
    cl_pair = []
    # search idx
    k = 0
    if cl_idx.shape[0]>0:
        for j in range(len(ends)-1):
            if (ends[j:j+2] == cl_idx[k]).all():
                cl_pair.append([j,j+1])
                k+=1
            if k==cl_idx.shape[0]:
                break

    if len(ends) > min_ends:
        cl_graph = {}
        for i in range(len(ends)):
            cl_graph[i] = set()
        for i in cl_pair:
            cl_graph[i[0]] = set([i[1],])
            cl_graph[i[1]] = set([i[0],])
        
        labels, mus = dpmeans(ends.reshape([-1,1]), w=w, cl_graph=cl_graph, Lambda=lbd, scores=scores).fit(ends.reshape([-1,1]))
    elif len(ends) == 0:
        continue
    else:
        mus = ends.mean()
        labels = np.zeros(ends.shape)

    k = len(mus.flatten())
    peaks.append(mus.flatten())
    pi[idx] = labels+n_peaks
    n_peaks += k

np.save(os.path.join(path, f'read_assignment'), pi)
np.save(os.path.join(path, f'called_nucs'), np.concatenate(peaks))
