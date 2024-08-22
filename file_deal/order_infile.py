#-*- encoding:utf-8 -*-
'''
将输入文件按照某列排队并输出前n行
'''

import argparse
import sys
import os
import pandas as pd 
bindir = os.path.abspath(os.path.dirname(__file__))

__author__='author'
__mail__= 'author@genome.cn'
def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument('-i', '--infile', help='infile', dest='infile', required = True )
    parser.add_argument('-o', '--outfile', help='outfile', dest='outfile', required = True )
    parser.add_argument('-col', '--column', help='column to sort', dest='column', required = True )
    parser.add_argument('-t', '--top', help='输出排序后的top多少，默认为0-全部输出 ', dest='top', required = True, type=int, default=0 )
    parser.add_argument('-a','--ascending',help='是否升序',dest='ascending',default=True)
    args=parser.parse_args()
    
    data = pd.read_csv(args.infile, sep="\t")
    data_sort = data.sort_values( by= args.column, ascending=args.ascending)
    if args.top == 0 :
        data_out = data_sort
    else:
        data_out = data_sort[0:args.top]
    data_out.to_csv( args.outfile, sep="\t", index=False)


if __name__ == '__main__':
    main()