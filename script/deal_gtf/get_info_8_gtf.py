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
def get_value(in_str, target ):
	if target in in_str:
		tmp = in_str.split('"')
		value = tmp[1]
		return value
	else:
		return 0
	

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	args=parser.parse_args()
	
	with open(args.input, 'r') as infile, open(args.output, 'w') as outfile :
		for line in infile:
			if line.startswith('#') : continue
			tmp = line.rstrip().split('\t')
			content = tmp[8]
			tt = content.split(';')
			gene_name = ''
			gene_id = ''
			for i in tt :
				if 'gene_id' in i :
					gene_id = get_value( i,'gene_id' )
					continue 
				if 'gene_name' in i :
					gene_name = get_value( i,'gene_name' )
			if gene_id and gene_name :
				outfile.write(gene_id+'\t'+gene_name+'\n')
			else:
				continue

if __name__ == '__main__':
	main()