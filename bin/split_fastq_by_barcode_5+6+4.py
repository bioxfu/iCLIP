#! /usr/bin/env python

import sys
import gzip

f = gzip.open(sys.argv[1], 'r')
output = gzip.open(sys.argv[2], 'w')

lst = sys.argv[2].split('.')
barcode = lst[1]
#print(barcode)

for line in f:
    fq_line1 = line
    
    # python2.x
    #fq_line2 = f.next()
    #fq_line3 = f.next()
    #fq_line4 = f.next()
    
    # python3.x
    fq_line2 = next(f)
    fq_line3 = next(f)
    fq_line4 = next(f)
    
    # python2.x
    # bc = fq_line2[5:11]
    
    # python3.x
    bc = fq_line2[5:11].decode("utf-8")
    
    # print(bc)
    if barcode == bc:
        #print('ok')
        output.write(fq_line1)
        output.write(fq_line2)
        output.write(fq_line3)
        output.write(fq_line4)
        
f.close()
output.close()

