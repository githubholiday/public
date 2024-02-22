
'''
功能：合并6个cluster的差异基因，如果在某个cluster中没有的，则fc=0
'''


import argparse
import sys # 获取命令行参数
import os
import re

bindir = os.path.dirname(os.path.abspath( __file__ ))

def store_gene( infile, gene_dict , cluster_list):
    with open(infile, 'r') as input:
        cluser_name = ''
        for line in input:
            tmp = line.rstrip().split('\t')
            if line.startswith('GeneName'):
                cluster_name = tmp[1]
                continue
            gene_name = tmp[0]
            fc = tmp[1]
            if gene_name not in gene_dict:
                gene_dict[gene_name] = {}
                for c in cluster_list :
                    gene_dict[gene_name][c] = '0'
            gene_dict[gene_name][cluster_name] = fc
    return gene_dict 
                
def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument( '-i', '--infile', help='de gene file of cmp1', dest='infile', required=True)
    parser.add_argument( '-o', '--outfile', help='outfile', dest='outfile', required=True)
    args = parser.parse_args()
    cluster_list = ['c0','c1','c2','c3','c4','c5','c6']
    gene_dict = {}
    for i in args.infile:
        gene_dict = store_gene(i, gene_dict, cluster_list )
    with  open( args.outfile, 'w') as output :
        head = ["#GeneName"]+cluster_list
        output.write('\t'.join(head)+'\n')
        for gene_name in gene_dict :
            gene_value = []
            for c in cluster_list :
                fc = gene_dict[gene_name][c]
                gene_value.append(fc)
            out_value = [gene_name]+gene_value
            output.write('\t'.join(out_value)+'\n')

if __name__ == '__main__':
    main()
