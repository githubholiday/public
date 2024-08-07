#! /usr/bin/env python3
import argparse
import sys
import os
import re
import time

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

def read_file( file ) :
    infile=open( file, 'r' )
    line =  infile.readlines()[1].split('\t')
    return line

def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument('-i','--input',help='input file',dest='input',nargs='+',required=True)
    parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
    args=parser.parse_args()

    title = 'Sample\tHost Genome Length(bp)\tTotal Reads\tTotal Bases(bp)\tMapped Host Reads\tMapped Host Rate(%)\n'
    args.output.write(title)
    for file in args.input :
        check_file_exists( file )
        sample = os.path.basename(file).split('.')[0]
        line = read_file( file )
        out = [ sample ] + line
        args.output.write('\t'.join(out) + '\n')


if __name__ == '__main__':
    main()
