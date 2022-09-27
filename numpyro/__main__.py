#!/usr/bin/env python

import gzip
import numpy as np
from sample import run_sampler

row_multiply = np.array([1, 2, 3], dtype=np.int8)

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def count_codons(file_name):
    with open(file_name, 'rb') as test_f:
        zipped = test_f.read(2) == b'\x1f\x8b'
    if zipped:
        fh = gzip.open(file_name, 'rt')
    else:
        fh = open(file_name, 'rt')
    with fh as fasta:
        for h, s in read_fasta(fasta):
                s = np.fromstring(s.lower(), dtype=np.int8)
                s[(s!=97) & (s!=99) & (s!=103) & (s!=116)] = 110
                codon_s = s.reshape(3, -1)
                # TODO Remove any Ns
                codon_s[s==97] = 0
                codon_s[s==99] = 1
                codon_s[s==103] = 2
                codon_s[s==110] = 3
                codon_idx_vec = np.sum(np.matmul(codon_s, row_multiply), 0)
                # TODO look this up in a table which maps to indicies in the X count



def main():
    # Need to read in data
    # see https://github.com/gtonkinhill/pairsnp-python/blob/master/pairsnp/pairsnp.py
    # and https://github.com/gtonkinhill/panaroo/blob/master/panaroo/prokka.py
    # but should be easy enough, just use a static lookup from codon to idx
    # (eventually reduce sum/merge would work if wanting to parallelise)
    print(run_sampler(X, pi_eq, warmup, samples))

if __name__ == "__main__":
    main()
