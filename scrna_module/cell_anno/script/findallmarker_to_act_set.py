#! /usr/bin/env python3
"""
	This script is used to convert the output of findallmarker
	输入： 
    p_val	avg_log2FC	pct.1	pct.2	p_val_adj	cluster	gene
    0.00e+00	7.115580	0.016	0.000	0.0004033	Naive CD4 T	GTSCR1
    0.00e+00	5.853999	0.022	0.001	0.0000122	Naive CD4 T	REG4
    0.00e+00	4.791871	0.034	0.003	0.0000001	Naive CD4 T	NOG
    8.71e-05	4.425698	0.010	0.001	1.0000000	Naive CD4 T	ST6GALNAC1
    8.67e-05	4.423395	0.010	0.001	1.0000000	Naive CD4 T	AXIN2

   输出： 
   cluster1:LTB,AQP3,CD2,IL32,IL7R,TNFRSF4,CD3D,LDHB,CORO1B,CD3E,TRAT1,TRADD,CD40LG,SPOCK2,OPTN,PPP2R5C,CD27,TTC39C,SUSD3,GIMAP4,AES,MAL,LAT,NOSIP,TNFAIP8,ARHGAP15,PRDX2,ANXA1,TNFAIP3,SLC2A3
    cluster2:S100A9,S100A8,LYZ,LGALS2,FCN1,TYROBP,CST3,MS4A6A,CD14,TYMP,LGALS1,GSTP1,GRN,CFD,AIF1,GPX1,LGALS3,AP1S2,CEBPD,CTSS,NPC2,CDA,NCF2,LST1,CD68,S100A6,FCGRT,CCL3,CFP,IL8
    cluster3:FCER1A,HLA-DPB1,HLA-DQA1,HLA-DRB1,HLA-DPA1,HLA-DRA,CST3,CLEC10A,HLA-DQB1,HLA-DRB5,ALDH2,CD74,HLA-DQA2,CD1C,CPVL,SERPINF1,LYZ,PLD4,HLA-DMA,CNDP2,ENHO,CFP,C1orf162,CD302,COMMD3,GDI2,AP1S2,GSTP1,CAPG,GRN

"""
import argparse
import sys
import os
import re
bindir = os.path.abspath(os.path.dirname(__file__))

__author__='Liu Tao'
__mail__= 'liutao@seqwisdom.com'

pat1=re.compile('^\s*$')

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',type=open,required=True)
	parser.add_argument('-n','--number',help='number of top gene',dest='number',type=int,default=20)
	parser.add_argument('-o','--output',help='output file',dest='output', required=True , type=argparse.FileType('w'))
	#parser.add_argument('-m','--mm',help='output file',dest='mm',action='store_false')
	args=parser.parse_args()

	headder_name_list = [ "cluster" , "gene" , "avg_log2FC" , "pct.1" , "pct.2" ,  "p_val_adj"]
	header_name_index = []
	
	result = {}

	for (count , line) in enumerate(args.input):
		if line.startswith('#') or re.search(pat1,line):continue
		tmp=line.rstrip().split('\t')
		if count == 0:
			for i in headder_name_list:
				a_index = tmp.index(i)
				header_name_index.append(a_index)
			print(header_name_index)
		else :
			cluster = tmp[header_name_index[0]]
			gene = tmp[header_name_index[1]]
			avg_log2FC = float(tmp[header_name_index[2]])
			pct1 = tmp[header_name_index[3]]
			pct2 = tmp[header_name_index[4]]
			pct_diff = float(pct1) - float(pct2)
			p_val_adj = float(tmp[header_name_index[5]])
			if p_val_adj > 0.05:continue
			if cluster not in result:
				result[cluster] = [[ gene , avg_log2FC , pct_diff  ]]
			else:
				result[cluster].append([ gene , avg_log2FC , pct_diff  ])
	
	for cluster in result:
		gene_info = result[cluster]
		sorted_gene_info1 = sorted(gene_info,key=lambda x:x[1],reverse=True)
		top_n_gene_name1 = [ x[0] for x in sorted_gene_info1[:args.number] ]
		print('{0}_avg_log2FC:{1}'.format(cluster,','.join(top_n_gene_name1)),file=args.output)
		
		sorted_gene_info2 = sorted(gene_info,key=lambda x:x[2],reverse=True)
		top_n_gene_name2 = [ x[0] for x in sorted_gene_info2[:args.number] ]
		print('{0}_pct_diff:{1}'.format(cluster,','.join(top_n_gene_name2)),file=args.output)


		top_n_gene_name12 = [ x[0] for x in sorted_gene_info1[:args.number*2] ]
		top_n_gene_name22 = [ x[0] for x in sorted_gene_info2[:args.number*2] ]
		top_n_gene_name_overlap = list(set(top_n_gene_name12).intersection(set(top_n_gene_name22)))
		print('{0}_overlap:{1}'.format(cluster,','.join(top_n_gene_name_overlap)),file=args.output)
			
			

			
			
			



			
			
			
			

			
		
	

if __name__ == '__main__':
	main()
