#! /usr/bin/env python

import sys

for line in sys.stdin:
    line = line.strip()
    ID = line[0:5]+line[11:15]
    barcode = line[5:11]
    seq = line[15:]
    #print barcode+'\t'+ID+'\t'+seq
    print(seq+'\t'+ID)
