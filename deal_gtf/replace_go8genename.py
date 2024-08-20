##### cat example.py 
#! /usr/bin/env python3
import argparse
import sys
import os
import re
import datetime
import pandas as pd

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
	parser.add_argument('-g','--goin',help='goin',dest='goin',required=True)
	parser.add_argument('-r','--relation',help='gene name and gene id relation file',dest='relation',required=True)
	parser.add_argument('-o','--gout',help='gout',dest='gout',required=True)
	args=parser.parse_args()
	relation_dict = get_relation(args.relation)
	with open(args.goin, 'r') as infile, open(args.gout, 'w') as outfile :
		for line in infile:
			tmp = line.rstrip().split('\t')
			gene_id = tmp[0]
			if gene_id in relation_dict :
				gene_name = relation_dict[gene_id]
				tmp[0] = gene_name
				outfile.write('\t'.join(tmp)+'\n')
			else:
				outfile.write('\t'.join(tmp)+'\n')


if __name__ == '__main__':
	main()