#! /usr/bin/env python

import sys

for line in sys.stdin:
    line = line.strip()
    ID = line[0:3]+line[7:9]
    barcode = line[3:7]
    seq = line[9:]
    #print barcode+'\t'+ID+'\t'+seq
    print(seq+'\t'+ID)
