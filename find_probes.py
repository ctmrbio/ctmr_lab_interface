#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from pathlib import Path
from sys import argv, exit
from collections import defaultdict


"""
For a rolling window on a given fasta sequence, count 
occurences of the same subsequence in other sequences

Useful for developing FISH/FACS probes
"""

__author__ = "CTMR, Luisa W. Hugerth"
__date__ = "2021"
__version__ = "0.1"



def find_substrings(target, length):
    """
    finds each substring of a given length in a sequence
    """
    sequence = target.seq
    L = len(sequence)

    allseqs = dict()

    if(L < length):
        return(sequence)
    else:
        for i in range(0, L-length):
            allseqs[str(i)] = str(sequence[i:i+length])

    return(allseqs)


def count_occurrences(all_records, target_id, candidates):
    """
    goes through each record and counts occurences of candidates
    """
    counts = defaultdict(int)
    for pos, probe in candidates.items():
        for ID, record in all_records.items():
            if ID != target_id:
                seq = record.seq
                if probe in seq:
                    counts[probe]+=1
    return(counts)


def print_output(candidates, occurrences):
    print("Position\tSequence\tOff-target matches")
    for pos, seq in candidates.items():
        if(seq in occurrences.keys()):
            print("\t".join([pos, seq, str(occurrences[seq])]))
        else:
            print("\t".join([pos, seq, str(0)]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', 
        required=True,
        help='Path to fasta file',
    )
    parser.add_argument('-t', '--target',
        required=True,
        help='ID of the target sequence'
    )
    parser.add_argument('-l', '--length',
        default=30,
        type=int,
        help='Length of the desired probe'
    )
    if len(argv) < 2:
        parser.print_help()
        exit()

    args = parser.parse_args()

    record_dict = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    candidates = find_substrings(record_dict[args.target], args.length)
    #print(candidates)
    occurrences = count_occurrences(record_dict, args.target, candidates)
    #print(occurrences)
    print_output(candidates, occurrences)
