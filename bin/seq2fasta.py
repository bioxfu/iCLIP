#! /usr/bin/env python

import sys

num = 0
for line in sys.stdin:
    l = line.strip()
    if len(l) >= 5:
        num += 1
        print('>'+str(num)+'\n'+line, end='')
    
