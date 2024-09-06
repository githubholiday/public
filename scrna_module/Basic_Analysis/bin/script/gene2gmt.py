
'''
将第一列为gene信息的基因列表转化为gene.gmt格式
'''

import argparse # 程序定义它需要的参数，从sys.argv解析出参数
import sys # 获取命令行参数
import os            #应用到系统
import time
import sys

bindir = os.path.dirname(os.path.abspath( __file__ ))

def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument( '-i', '--infile', help='dir to check', dest = 'infile', nargs="+")
    #parser.add_argument( '--header', help='header in file or not', dest = 'header')
    parser.add_argument( '-o','--outfile', help='header in file or not', dest = 'outfile')
    args = parser.parse_args()

    for infile in args.infile:
        infile_abspath = os.path.abspath(infile)
        gene_set_name = infile_abspath.split('/')[-1].split('.')[0]
        with open(infile, 'r') as input, open( args.output, 'w') as output:
            gene_set_info = [gene_set_name,gene_set_name]
            for line in input:
                tmp = line.rstrip().split('\t')
                gene_id = tmp[0]
                gene_set_info.append(gene_id)
        output.write('\t'.join(gene_set_info)+'\n')


if __name__ == '__main__':
    main()
