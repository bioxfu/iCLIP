#! /usr/bin/env python

import sys
import os

if len(sys.argv) < 2:
    print "Usage: fastqc_html_reads_stat.py fastqc.html"
else:
    cmd = r"grep 'Filename' %s|sed -r 's/.+Filename<\/td><td>//'|sed -r 's/<\/td><\/tr><tr><td>Sequences flagged.+//'|sed -r 's/<.+>/\t/'" % sys.argv[1] 
    os.system(cmd)

