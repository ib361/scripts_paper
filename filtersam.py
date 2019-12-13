#!/bin/usr/env python 3
import sys

if len(sys.argv[1:])<2:
    sys.stderr.write('\nERROR: missing argument\n\n')
    sys.stderr.write('Usage: script.py /path/to/*mitioonly.sam <stats_filename.tsv> <readID_list_filename.tsv>\n\n')
    sys.exit(1)

path, stats_outfile, readids_outfile = sys.argv[1:]

print('Importing modules...\n')
import fileinput
from glob import glob
import numpy as np
import pandas as pd

a = sorted(glob(path+'*mitoonly.sam'))

def parse_sa_flag(saflag):
    saflag = saflag.split(',')
    chrom = saflag[0].split(':')[-1]
    strand = saflag[2]
    v = [chrom, strand]
    return v

def check_read_origin(i,chrmt_name):
    l = len(i)
    newi = np.unique(np.array(i),axis=0)
    s = sum(newi[:,0]==chrmt_name)
    if s == len(newi):
        strand = newi[:,1].tolist()
        if '+' in strand and '-' in strand:
            tag = 'MITO_CHIMERIC'
        else:
            if l > 2 :
                tag = 'MITO_CHIMERIC' #cases when more than two alignments of the same read happen to be on MT
            else:
                tag = 'MITO'
    else:
        tag = 'NUMT'
    return tag

df_list = []
toexcludelist = []
print('Calculating stats on mitoonly sam files...\n')
for sam in a:
    dic = {}
    readcount = 0
    for line in fileinput.input(sam):
        line = line.split("\t")
        if len(line) == 20:
            rid = line[0]
            readcount += 1
        if len(line) ==21:
            rid = line[0]
            saflag = line[20].strip()
            saflag = saflag.split(';')
            saflag = list(filter(None,saflag)) #remove empty values from list
            if rid not in dic:
                dic[rid] = list(map(lambda x:parse_sa_flag(x),saflag))
                readcount += 1
            else:
                dic[rid].extend(list(map(lambda x:parse_sa_flag(x),saflag)))
    fileinput.close()
    #populate a new dictionary with tags defining the origin of each read
    dict_tags = {}
    for i in dic:
        dict_tags[i] = check_read_origin(dic[i],'MT')
    #print stats
    mito = list(dict_tags.values()).count('MITO')
    numt = list(dict_tags.values()).count('NUMT')
    mito_chimeric = list(dict_tags.values()).count('MITO_CHIMERIC')
    unique = readcount - (mito + numt + mito_chimeric)
    tmp = pd.DataFrame([['MITO (non SA)',unique],['MITO', mito], ['NUMT',numt], ['MITO CHIMERIC', mito_chimeric]])
    tmp.columns = ['Tag','Counts']
    tmp.insert(2,'Perc',(tmp['Counts']/readcount)*100)
    df_list.append(tmp)
    numtlist = [key for (key,value) in dict_tags.items() if value == 'NUMT']
    chimericlist = [key for (key,value) in dict_tags.items() if value == 'MITO_CHIMERIC']    
    toexcludelist.append(numtlist)
    toexcludelist.append(chimericlist)

flat_toexclude = []
for sublist in toexcludelist:
    for item in sublist:
        flat_toexclude.append(item)

df_final = pd.concat(df_list)
print('Saving results into file...\n')
df_final.to_csv(stats_outfile, sep = '\t', index = None)
flat_toexclude_df = pd.DataFrame(flat_toexclude)
flat_toexclude_df.to_csv(readids_outfile, sep = '\t', index = None, header = None)
