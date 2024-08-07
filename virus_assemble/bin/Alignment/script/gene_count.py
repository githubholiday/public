'''
将每个基因比对上的reads进行统计
'''
#! /usr/bin/env python2.7
import argparse
import sys
import os
import re
import subprocess

__author__='Liu Tao'
__mail__= 'taoliu@annoroad.com'

pat1=re.compile('^\s+$')

def defineRead(flag):
	bin_count = bin(flag).replace("0b",'')
	try:
		if bin_count[-3] == '1':
			return '3'
		elif bin_count[-7] == '1' :
			return '1'
		elif bin_count[-8] == '1':
			return '2'
	except:
		pass

def read_bam(bam, samtools):
	r_dict={}
	proc = subprocess.Popen([samtools,'view',bam],stdout = subprocess.PIPE)
	for line in proc.stdout:
		if line.startswith(b'#'): continue
		tmp=line.rstrip().decode().split("\t")
		name, flag , chr, qual = tmp[0], int(tmp[1]), tmp[2], tmp[4]
		read_number = defineRead(flag)
		if read_number == '3' : continue
		if not chr in r_dict:
			r_dict[chr]=0
		r_dict[chr]+=1
	return r_dict

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file,bam format',dest='input',required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	parser.add_argument('-s','--samtools',help='samtools path',dest='samtools',required=True)
	args=parser.parse_args()

	d_raw_count=read_bam(args.input, args.samtools)
	for contig in sorted(d_raw_count):
		count = d_raw_count[contig]
		args.output.write('{0}\t{1}\n'.format(contig,count))

if __name__=='__main__':
	main()
