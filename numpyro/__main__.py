#!/usr/bin/env python

from sample import run_sampler

def main():
    # Need to read in data
    # see https://github.com/gtonkinhill/pairsnp-python/blob/master/pairsnp/pairsnp.py
    # and https://github.com/gtonkinhill/panaroo/blob/master/panaroo/prokka.py
    # but should be easy enough, just use a static lookup from codon to idx
    # (eventually reduce sum/merge would work if wanting to parallelise)
    print(run_sampler(X, pi_eq, warmup, samples))

if __name__ == "__main__":
    main()