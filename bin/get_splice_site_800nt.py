#! /usr/bin/env python

import sys
import re

intronLen = 800 - 40

for line in sys.stdin:
    if not line.startswith('\t'):
        line_lst = line.strip().split('\t')
        temp = line_lst[0].split('|')
        if len(temp) == 10:  ## old version of AS output
            chrom, s1, s2, s3, s4, s5, s6, strand, name = temp[2],temp[4],temp[5],temp[6],temp[7],temp[8],temp[9],temp[3],temp[1]
        elif len(temp) == 7: ## new version of AS output
            chrom, s2, s3, s4, s5, strand, name = temp
            s1 = int(s2) - 100
            s6 = int(s5) + 100
        elif len(temp) == 8: ## AS output from rMAT
            chrom, strand = temp[:2]
            s1, s2, s3, s4, s5, s6 = sorted([int(x) for x in temp[2:]])
                            
        s1, s2, s3, s4, s5, s6= [int(x) for x in [s1, s2, s3, s4, s5, s6]]
        change = line_lst[1]
        region1 = [s2-40, s2+intronLen]
        region2 = [s3-intronLen,s3+40]
        region3 = [s4-40,s4+intronLen]
        region4 = [s5-intronLen,s5+40]
        if s2-s1 < 40*2:
            region1[0] = int(s2-(s2-s1)/2)
        if s3-s2 < intronLen*2:
            region1[1] = int(s2+(s3-s2)/2)
            region2[0] = int(s3-(s3-s2)/2)
        if s4-s3 < 40*2:
            region2[1] = int(s3+(s4-s3)/2)
            region3[0] = int(s4-(s4-s3)/2)
        if s5-s4 < intronLen*2:
            region3[1] = int(s4+(s5-s4)/2)
            region4[0] = int(s5-(s5-s4)/2)
        if s6-s5 < 40*2:
            region4[1] = int(s5+(s6-s5)/2)
        
        for i, j in enumerate(range(s2,region1[0],-10)):
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j-10, j, line_lst[0], change, strand, 1, -i))
        for i, j in enumerate(range(s2,region1[1],10)):
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j, j+10, line_lst[0], change, strand, 1, i+1))

        for i, j in enumerate(range(s3,region2[0],-10)):
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j-10, j, line_lst[0], change, strand, 2, -i))
        for i, j in enumerate(range(s3,region2[1],10)):
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j, j+10, line_lst[0], change, strand, 2, i+1))
            
        for i, j in enumerate(range(s4,region3[0],-10)):
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j-10, j, line_lst[0], change, strand, 3, -i))
        for i, j in enumerate(range(s4,region3[1],10)):
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j, j+10, line_lst[0], change, strand, 3, i+1))
               
        for i, j in enumerate(range(s5,region4[0],-10)):
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j-10, j, line_lst[0], change, strand, 4, -i))
        for i, j in enumerate(range(s5,region4[1],10)):
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, j, j+10, line_lst[0], change, strand, 4, i+1))
            
