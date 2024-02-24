#!/usr/bin/env python

import re, string, itertools, os
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO, bgzf
from io import StringIO
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("read_input")
parser.add_argument("read_output")

args = parser.parse_args()

## The trimming function
def trim_hps(title, seq, qual):
    ## split consecutive elements
    rep_list = re.findall(r'((.)\2*)', seq)
    # Get first and last repetitive sets
    first_rep = rep_list[0][0]
    last_rep = rep_list[-1][0]
    # check if there are more than 1 or 2 repeat sets
    if (len(first_rep) > 1 or len(last_rep) > 1) and len(rep_list) > 2:
        ## trim the sets
        if len(first_rep) > 1:
            first_rep = rep_list[0][1]

        if len(last_rep) > 1:
            last_rep = rep_list[-1][1]

        middle_portion = ''.join([x[0] for x in rep_list][1:-1])

        new_seq = ''.join(map(str, (first_rep, middle_portion, last_rep)))

        # adjust the quality to the same length
        qual_start = len(rep_list[0][0])-len(rep_list[0][1])
        qual_stop = len(rep_list[-1][0])-len(rep_list[-1][1])
        new_qual = qual[qual_start: len(qual)-qual_stop]
        #save the record
        fastq_string = "@%s\n%s\n+\n%s\n" % (title, new_seq, new_qual)
        record = SeqIO.read(StringIO(fastq_string), "fastq")

    else:
        fastq_string = "@%s\n%s\n+\n%s\n" % (title, seq, qual)
        record = SeqIO.read(StringIO(fastq_string), "fastq")

    return record

## applied looping through reads
count = 0
with gzip.open(args.read_input, "rt") as in_handle:
    with bgzf.BgzfWriter(args.read_output, "wb") as out_handle:
        for title, seq, qual in FastqGeneralIterator(in_handle):
            SeqIO.write(trim_hps(title, seq, qual), out_handle, "fastq")
            count += 1


print("Saved %i records from %s to %s" % (count, args.read_input, args.read_output))
