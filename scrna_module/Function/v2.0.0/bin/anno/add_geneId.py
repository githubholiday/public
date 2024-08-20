'''

升级：
2024-08-13 tcf
1. 读取的-i参数中的第一列Gene列中的基因名根据-分割，如果最后一个字符是数字，则直接去掉，如果不是字符忽略
原因：从rds中获取的基因名会增加-1-2，导致其找不到对应的gene_id，而被过滤掉
'''
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
			if gene_name not in relation_dict:
				relation_dict[gene_name] = gene_id
	return relation_dict

	

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--in',help='infile',dest='infile',required=True)
	parser.add_argument('-r','--relation',help='gene name and gene id relation file',dest='relation',required=True)
	parser.add_argument('-o','--gout',help='gout',dest='gout',required=True)
	args=parser.parse_args()
	relation_dict = get_relation(args.relation)
	with open(args.infile, 'r') as infile, open(args.gout, 'w') as outfile :
		for line in infile:
			if line.startswith("#") or line.startswith('Gene') : 
				tt = line.replace('Gene',"GeneName")
				outfile.write("Gene\t"+tt)
			tmp = line.rstrip().split('\t')
			gene_name = tmp[0]
			if gene_name in relation_dict :
				gene_id = relation_dict[gene_name]
				tmp = [gene_id] + tmp
				outfile.write('\t'.join(tmp)+'\n')
			else:
				#如果基因名匹配不上的话，则将基因名称后面的-1-2去掉
				gene_name_t = tmp[0]
				gene_name_list = gene_name_t.split("-")
				if gene_name_list[-1].isdigit():#判断最后一个元素是否为数字
					gene_name_list.pop()
					gene_name = '-'.join(gene_name_list) #将最后一个元素去掉
				else:
					gene_name = gene_name_t
				if gene_name in relation_dict :
					gene_id = relation_dict[gene_name]
					tmp = [gene_id] + tmp
					outfile.write('\t'.join(tmp)+'\n')
				else:
				#tmp = ['']+tm		p
					print("{0} has no relative gene id".format( gene_name_t ))
				#outfile.write('\t'.join(tmp)+'\n')


if __name__ == '__main__':
	main()
