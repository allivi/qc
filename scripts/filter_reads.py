#!/usr/bin/env python
#USAGE: python filter_reads.py [include|exclude] reads.fq filter_ids.py
import sys

include = True if sys.argv[1] == 'include' else False

ids_f = open(sys.argv[3])
ids = [l.strip() for l in ids_f]
reads_f = open(sys.argv[2])
reads = []

reads_iter = iter(reads_f)
for l in reads_iter:
    read_id = l[1:].strip()
    need_print = False
    if read_id in ids:
        need_print = include
    else:
        need_print = not include
    if need_print:
        print l.strip()
        for i in range(3):
            try:
                print (next(reads_iter).strip())
            except StopIteration:
                sys.exit()
