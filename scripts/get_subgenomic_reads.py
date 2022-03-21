import os
import sys
import numpy as np
import pysam
import pandas as pd
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument('-b','--bam',dest='bam',
                    help="""input bam file(s) """)
parser.add_argument('-c','--chrom',dest='chrom',
                    help="""relevant chromosome""")
parser.add_argument('-o','--out',dest='out',
                    help="""output csv file""")

args=parser.parse_args()

chroms={'SCoV2': 'NC_045512.2',
        'SCoV1': 'NC_004718.3',
        'MERS': 'NC_019843.3',
        'NL63': 'JX504050.1'}

if args.chrom in chroms:
    chrom=chroms[args.chrom]

n_subgenomic={}
n_viral={}
n_tot={}
for bamfile in args.bam.split(','):
    print('reading '+bamfile)
    n_subgenomic[bamfile]=0
    n_viral[bamfile]=0
    with pysam.Samfile(bamfile) as bam:
        n_tot[bamfile]=bam.mapped
        if chrom in bam.references:
            for read in bam.fetch(chrom):
                n_viral[bamfile]+=1
                # look for reads that originate in the first 120 nt and have exactly one gap of at least 10nt at least 15 nt away from the ends
                dp=np.diff(read.positions)
                gaps=dp[dp > 1]
                gappos=np.where(dp > 1)[0]
                if read.reference_start < 120 and len(gaps)==1 and gaps[0] > 10 and gappos[0] > 15 and gappos[0] < read.qlen-15:
                    n_subgenomic[bamfile]+=1
                    pos=read.positions
                    dpp=gappos
                    esize=[pos[dpp[0]]-pos[0]+1]+[pos[dpp[i]]-pos[dpp[i-1]+1]+1 for i in range(1,len(dpp))] + [pos[-1]-pos[dpp[-1]+1]+1]
                    estart=[0]+[pos[p+1]-pos[0] for p in dpp[:-1]] + [pos[dpp[-1]+1]-pos[0]]
                    if pos[0]+estart[-1]+esize[-1]!=pos[-1]+1:
                        raise Exception('stop')
                    #print('{0}\t{1}\t{2}\t{3}\t0\t.\t{1}\t{2}\t{7}\t{4}\t{5}\t{6}'.format(chrom,pos[0],pos[-1]+1,read.qname,len(esize),','.join(map(str,esize))+',',','.join(map(str,estart))+',','0,0,0'))

pd.DataFrame.from_dict({'n_tot': n_tot, 'n_viral': n_viral, 'n_subgenomic': n_subgenomic},orient='columns').to_csv(args.out)
