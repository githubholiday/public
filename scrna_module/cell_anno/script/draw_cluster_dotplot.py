#! /usr/bin/env python3
import argparse
import sys
import os
import re
import scanpy as sc
import matplotlib.pyplot as plt

bindir = os.path.abspath(os.path.dirname(__file__))

__author__='Liu Tao'
__mail__= 'liutao@seqwisdom.com'

pat1=re.compile('^\s*$')


def parse_gmt(gmt):
	gmt_dict={}
	for line in gmt:
		if line.startswith('#') or re.search(pat1,line):continue
		tmp=line.rstrip().split('\t')
		gmt_dict[tmp[0]]=tmp[2:]
	return gmt_dict

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input rds',dest='input', required=True)
	parser.add_argument('-g','--gmt',help='gmt file',dest='gmt',type=argparse.FileType('r'),required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	parser.add_argument('-c','--cluster',help='cluster column',dest='cluster',default="harmony_clusters")
	#parser.add_argument('-m','--mm',help='output file',dest='mm',action='store_false')
	args=parser.parse_args()
    
	gmt_dict = parse_gmt(args.gmt)
	adata = sc.read_h5ad(args.input)
	sc.pl.dotplot(adata,  gmt_dict, dendrogram=True , groupby= args.cluster)
	plt.savefig(args.output)
	
	

if __name__ == '__main__':
	main()
