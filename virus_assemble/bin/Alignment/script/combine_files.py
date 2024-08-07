'''
function :将两个文件（具有相似格式的,至少有一列是相互包含的）,将file1和file2中指定列中的相同内容的行进行合，默认通过两个文件的第一列进行行合并；
参数说明:
	-i1 : 输入的第一个文件；（输出时会先输出文件1）
	-i2 : 输入的第二个文件；
	-c1 : 想要通过file1的哪一列进行数据合并，为整数，1,2,3等,从1开始计数,默认为第一列
	-c2 : 想要通过file2的哪一列进行数据合并，为整数，1,2,3等,从1开始计数，默认为第一列
	-o : 输出文件;输出的第一列数据是进行数据合并的列，后面依次是file1,file2的其余列
示例：
python3 combine_2files.py -i1 file1.txt -i2 file2.txt -c1 2 -c2 1 -o file1_2.txt
将file1和file2按照文件的第二列内容进行合并
'''
#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import sys
import os
import re
import time
bindir = os.path.abspath(os.path.dirname(__file__))
filename = os.path.basename(__file__)

__author__='tu chengfang'
__mail__= 'chengfangtu@annoroad.com'

pat1=re.compile('^\s+$')

def std( level, message ):
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    string = '{0} - {1} - {2} - {3}\n'.format( now, filename, level, message )
    if level == 'ERROR':
        sys.stderr.write( string )
    else:
        sys.stdout.write( string )
	
def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i1','--input1',help='input file1',required=True)
	parser.add_argument('-i2','--input2',help='input file2',type=open,required=True)
	parser.add_argument('-c1','--col1',help='the col you want to combine by',type=int,default=1)
	parser.add_argument('-c2','--col2',help='the col you want to combine by',type=int,default=1)
	parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	args=parser.parse_args()
	
####----transfer the input1 to dic------####
	file1_dic = {}
	col1 = args.col1 - 1
	col2 = args.col2 - 1
	file1_target_list = []
	with open( args.input1 , 'r') as infile1 :
		for line in infile1 :
			if line.startswith('#') : continue 
			tmp = line.rstrip().split('\t')
			if len( tmp ) < col1 : continue
			file1_target = tmp[ col1 ]
			file1_target_list.append( file1_target )
			del tmp[ col1 ]
			tt = []
			for i in tmp :
				if not i :
					tt.append('**')
				else :
					tt.append( i )
			file1_dic[ file1_target ] = tt

	find_target = []
	for line in args.input2 :
		if re.search(pat1,line) : continue 
		lines = line.rstrip().split('\t')
		target = lines[col2]
		del lines[col2]
		if target in file1_target_list :
			if target in find_target :
				std("Warinning",'File2 already has {0} '.format( target ))
			else :
				find_target.append( target )
			file1_value = file1_dic[target]
			args.output.write( target+ '\t'+ '\t'.join(file1_value)+ '\t'+ '\t'.join( lines )+'\n' )
	
	find_target_set= set(find_target)
	missInfile1 = list(set(file1_target_list) - set(find_target_set))
	if len( missInfile1 ) == 0 :
		std("INFO",'the target of file1 all find in file2')
	else :
		std("Warinning",'the target of file1 is not find in file2 :{0}'.format( ','.join(list(missInfile1))))
if __name__ == '__main__':
	main()