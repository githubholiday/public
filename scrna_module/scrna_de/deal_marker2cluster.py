
'''
功能：
'''


import argparse
import sys # 获取命令行参数
import os
import re

bindir = os.path.dirname(os.path.abspath( __file__ ))

def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument( '-i', '--infile', help='de gene file of cmp1', dest='infile', required=True)
    #parser.add_argument( '-i2', '--infile2', help='de gene file of cmp2', dest='infile2', required=True)
    parser.add_argument( '-o', '--outfile', help='outfile', dest='outfile', required=True)
    parser.add_argument( '-p', '--prefix', help='outdir of result', dest='prefix', required=True)
    args = parser.parse_args()
    with open(args.infile, 'r') as input, open( args.outfile, 'w') as output :
        for line in input:
            if line.startswith(','):
                output.write("GeneName\t{0}\n".format(args.prefix))
                continue
            tmp = line.rstrip().split(',')
            gene_name = tmp[0]
            fc = tmp[2]
            output.write('\t'.join([gene_name, fc ])+'\n')





if __name__ == '__main__':
    main()
