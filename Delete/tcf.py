#-*- encoding:utf-8 -*-

import argparse
import sys
import os
import time
from datetime import datetime
import glob
bindir = os.path.abspath(os.path.dirname(__file__))


__author__='chengfangtu'
__mail__= 'chengfangtu@genome.cn'

def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument('-i','--infile',help='infile',dest='infile')
    parser.add_argument('-o','--outfile',help='outfile',dest='outfile')
    parser.add_argument('-o2','--outfile2',help='outfile2',dest='outfile2')
    args=parser.parse_args()
    import glob
    with open(args.infile, 'r') as input, open(args.outfile, 'w') as output,open( args.outfile2, 'w') as output2:
        for line in input:
            tmp = line.rstrip().split('\t')
            filter_path = tmp[4]
            filter_a_path = '{0}/Analysis'.format(filter_path)
            if os.path.exists(filter_a_path):
                output2.write('\t'.join(tmp)+'\n')
                continue
            filter_a_path_glob = glob.glob(filter_a_path+'/*')
            size_list = []
            for each_path in filter_a_path_glob:
                if each_path.endswith("report"):continue
                if each_path.endswith("STAT_result.xls"):continue
                if each_path.endswith("STAT_result.xls.mysql"):continue
                path_size = os.path.getsize(each_path)
                print(path_size)
                if path_size:
                    size_list.append('1')
            if size_list :
                output.write('\t'.join(tmp)+'\n')
            else:
                output2.write('\t'.join(tmp)+'\n')



if __name__ == '__main__':
    main()


