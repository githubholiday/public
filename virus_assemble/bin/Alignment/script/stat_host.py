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
    '''
    从report中获取reads数，单端reads数
    '''
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
    parser.add_argument('-t','--type',help='Mapped or Unmapped',dest='type',required=True)
    args=parser.parse_args()
    #如果输入的report是未必对上的，则-t 给定 unmap, 如果为比对上的，则为 map

    if args.type not in ["map","unmap"]:
        print("-t Please input map|unmap")
        sys.exit(1)
    total_reads = get_total( args.input )
    #从report里面获取的是unmaped
    report_reads = get_unmap( args.report) #int
    other_reads = int(total_reads) - report_reads
    if args.type == "map":
        mapped_reads = report_reads
        unmapped_reads = other_reads
    else:
        unmapped_reads = report_reads
        mapped_reads = other_reads
    
    map_rate = '{0:.2f}'.format( 100*mapped_reads/int(total_reads))
    unmap_rate = '{0:.2f}'.format( 100*unmapped_reads/int(total_reads))

    args.out.write("Sample\t"+args.sample+"\n")
    args.out.write("Total Reads\t"+total_reads+"\n")
    args.out.write("Mapped Reads\t"+str(mapped_reads)+"\n")
    args.out.write("Mapped Rate(%)\t"+map_rate+"\n")
    args.out.write("Unmapped Reads\t"+str(unmapped_reads))
    args.out.write("Unmapped Rate(%)\t"+str(unmap_rate))

if __name__ == '__main__':
    main()
