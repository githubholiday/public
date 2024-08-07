#! /usr/bin/env python3
import argparse
import sys
import os
import re
import time
import pandas as pd

__author__='zhang yue'
__mail__= 'yuezhang@genome.cn'

pat1=re.compile('^\s+$')
LOG = os.path.basename(__file__)

def std( level, message ):
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    string = '{0} - {1} - {2} - {3}\n'.format( now, LOG, level, message )
    if level == 'ERROR':
        sys.stderr.write( string )
    else:
        sys.stdout.write( string )

def check_file_exists( *file_list ) :
    for file in file_list :
        if os.path.exists( file ) :
            std( 'INFO', 'file : {0}'.format( file ) )
        else :
            std( 'ERROR', 'file is not exists : {0}'.format( file ) )

def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument('-i','--input',help='clean.bam file, removed host',dest='input',required=True)
    parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
    args=parser.parse_args()

    check_file_exists( args.input )
    df = pd.read_csv(args.input, low_memory=False, index_col='Tax', sep="\t")
    for c in df.columns :
        df[c] = df[c]/df[c].sum() * 100
    df.to_csv(args.output,sep='\t',index=True,header=True)


if __name__ == '__main__':
    main()
