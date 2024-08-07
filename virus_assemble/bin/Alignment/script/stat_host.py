#! /usr/bin/env python3
import argparse
import sys
import os
import re
import time

__author__='zhang yue'
__mail__= 'yuezhang@genome.cn'
__update__ = "liaorui"

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

def host_length( host ) :
    length = 0
    for line in host :
        if line.startswith('>') : continue
        length += len( line.rstrip() )
    return length

def get_total( stat ):
    with open(stat,'r') as infile:
         for line in infile:
             if '+ 0 primary' in line:
                total = line.split("+")[0].strip()
                return total


def get_unmap( report ):
    with open( report, 'r') as infile:
        for line in infile:
            if line.startswith('#total_reads'):
                unmaped_reads = line.split('\t')[1]
                return int(unmaped_reads)
            else: continue


def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument('-i','--input',help='clean.bam file, removed host',dest='input',required=True)
    parser.add_argument('-r','--report',help='*.report',dest='report',required=False)
    parser.add_argument('-o','--output',help='output file',dest='out',type=argparse.FileType('w'),required=True)
    parser.add_argument('-s','--sample',help='sample file',dest='sample',required=True)
    args=parser.parse_args()


    total_reads = get_total( args.input )
    unmaped_reads = get_unmap( args.report) #int
    mapped_reads = int(total_reads) - unmaped_reads
    map_rate = '{0:.2f}'.format( 100*mapped_reads/int(total_reads))

    args.out.write("Sample\t"+args.sample+"\n")
    args.out.write("Total Reads\t"+total_reads+"\n")
    args.out.write("Mapped Host Reads\t"+str(mapped_reads)+"\n")
    args.out.write("Mapped Host Rate(%)\t"+map_rate+"\n")
    args.out.write("Unmapped Host Reads\t"+str(unmaped_reads))

if __name__ == '__main__':
    main()
