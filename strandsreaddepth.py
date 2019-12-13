#!/bin/usr/env python3

import sys

if len(sys.argv[1:])<1:
    sys.stderr.write('\nERROR: missing argument\n\n')
    sys.stderr.write('Usage: script.py /path/to/nosa.sam/files <outfile_name.txt>\n\n')
    sys.exit(1)

path, outfile = sys.argv[1:]

print('Importing libraries...\n')
import pandas as pd
import glob
import fileinput

print('Counting strands...\n')
a = sorted(glob.glob(path+"*.sam"))
df_list = []
for i in a:
    dict_strand = {0:0, 16:0}
    for line in fileinput.input(i):
        if line.startswith("@"):
            pass
        else:
            line = line.split("\t")
            if line[2] == "MT":
                flag = int(line[1])
                dict_strand[flag]+=1
    plus = dict_strand[0]
    minus = dict_strand[16]
    total = plus + minus
    tmp = pd.DataFrame([['Plus', plus], ['Minus', minus], ['Total_reads', total]])
    tmp = tmp.set_index(0).transpose()
    tmp.insert(3, 'Perc_plus', round((plus/total)*100, 2))
    tmp.insert(4, 'Perc_minus', round((minus/total)*100, 2))
    df_list.append(tmp)
fileinput.close()
df_final = pd.concat(df_list)
print('Saving results into file...\n')
df_final.to_csv(outfile, sep = '\t', index = None, header = True)
