#! /usr/bin/python

import sys
bed_file = sys.argv[1]
bin_num = int(sys.argv[2])

for line in open(bed_file):
    line_lst = line.strip().split('\t')
    chrom, start, end, strand = line_lst[0], int(line_lst[1]), int(line_lst[2]), line_lst[5]
    regionLen = end - start + 1
    binLen = regionLen / bin_num
    bin_start = start
    bin_end = start
    i = 0
    while bin_end < end-regionLen%bin_num:
        i += 1
        bin_end = bin_start + binLen
        if strand == "+":
            print '%s\t%s\t%s\t%s\t%s' % (chrom, bin_start, bin_end, "\t".join(line_lst[3:]), i)
        else:
            print '%s\t%s\t%s\t%s\t%s' % (chrom, bin_start, bin_end, "\t".join(line_lst[3:]), bin_num+1-i)
        bin_start = bin_end
        
