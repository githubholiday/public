##### cat example.py 
#! /usr/bin/env python3
import argparse
import sys
import os
import re
import datetime

bindir = os.path.abspath(os.path.dirname(__file__))
filename=os.path.basename(__file__)

__author__='tu chengfang '
__mail__= 'chengfangtu@genome.cn'


class Log():
	def __init__( self, filename, funcname = '' ):
		self.filename = filename 
		self.funcname = funcname
	def format( self, level, message ) :
		date_now = datetime.datetime.now().strftime('%Y%m%d %H:%M:%S')
		formatter = ''
		if self.funcname == '' :
			formatter = '\n{0} - {1} - {2} - {3} \n'.format( date_now, self.filename, level, message )
		else :
			
			formatter = '\n{0} - {1} - {2} -  {3} - {4}\n'.format( date_now, self.filename, self.funcname, level, message )
		return formatter
	def info( self, message ):
		formatter = self.format( 'INFO', message )
		sys.stdout.write( formatter )
	def debug( self, message ) :
		formatter = self.format( 'DEBUG', message )
		sys.stdout.write( formatter )
	def warning( self, message ) :
		formatter = self.format( 'WARNING', message )
		sys.stdout.write( formatter )
	def error( self, message ) :
		formatter = self.format( 'ERROR', message )
		sys.stderr.write( formatter )
	def critical( self, message ) :
		formatter = self.format( 'CRITICAL', message )
		sys.stderr.write( formatter )
def get_relation( relation_file ):
	relation_dict = {}
	with open(relation_file, 'r') as infile:
		for line in infile:
			tmp = line.rstrip().split('\t')
			gene_id = tmp[0]
			gene_name = tmp[1]
			if gene_id not in relation_dict:
				relation_dict[gene_id] = gene_name
	return relation_dict

	

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--in',help='infile',dest='infile',required=True)
	parser.add_argument('-r','--relation',help='gene name and gene id relation file',dest='relation',required=True)
	parser.add_argument('-cn','--colname',help='colname of gene',dest='colname',required=True)
	parser.add_argument('-o','--out',help='out',dest='out',required=True)
	args=parser.parse_args()
	
	relation_dict = get_relation(args.relation) #{gene_id:gene_name}
	col_index = 0 
	print(args.colname)
	with open(args.infile, 'r') as infile, open(args.out, 'w') as outfile :
		for line_index, line in enumerate(infile):
			if line_index == 0 :
				tmp = line.rstrip().split('\t')
				for i_index, i in enumerate(tmp) :
					print(i)
					if i == args.colname:
						col_index = i_index
						print(col_index)
						outfile.write('\t'.join(tmp)+'\tGeneName\n')
						continue
				continue
			tmp = line.rstrip().split('\t')
			gene_id_str = tmp[col_index]
			gene_id_list = gene_id_str.split('/')
			gene_name_list = []
			for gene_id in gene_id_list:
				if gene_id not in relation_dict:
					#print("{0} not in relation_dict".format( gene_id ))
					continue
				gene_name = relation_dict[gene_id]
				gene_name_list.append(gene_name)
			gene_name_str = '/'.join(gene_name_list)
			outfile.write('\t'.join( tmp ) +'\t'+gene_name_str+'\n')

if __name__ == '__main__':
	main()

