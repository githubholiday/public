#!/usr/bin/env python3
import os
import sys 
import re
import argparse
bindir = os.path.abspath(os.path.dirname(__file__))

__author__='Yang Zhang'
__mail__= 'yangzhang@genome.cn'
__doc__='this file is used to append TPM from file2  for specified col in file1.' 


pat1=re.compile('^\s+$')

def read_annotation(f_file,col,tag ,r_dict):
	header = ''
	for count,  line in enumerate(f_file):
		if line.startswith('#') or re.search(pat1,line):continue
		tmp=line.rstrip().split('\t')
		id = tmp[col]
		ar_len=len(tmp)-1
		if count == 0 :
			if tag == 'ref':
				header = '\t'.join(line.rstrip().split('\t')[1:])
			elif tag == 'novel':
				header = "\t".join(tmp[1:])
		else:
			if not id in r_dict:
				r_dict[id] = ''
			if tag == 'ref':
				r_dict[id] = '\t'.join(line.rstrip().split('\t')[1:])
			elif tag == 'novel':
				r_dict[id] = "\t".join(tmp[1:])
	return r_dict,header,ar_len

def append_annotation(f_file,col,annotation,o_file,symbols):
	for line in f_file:
		if line.startswith('#') or re.search(pat1,line):
			o_file.write('{0}\t{1}\n'.format(line.rstrip(),'annotation'))
			continue
		if line.startswith('ID'):
			o_file.write('{0}\t{1}\n'.format(line.rstrip(),'annotation'))
			continue
		tmp=line.rstrip().split('\t')
		id = tmp[col]
		anno = [symbols]
		if id in annotation:
			anno = annotation[id]
		o_file.write('{0}\t{1}\n'.format(line.rstrip(),'\t'.join(anno[1:])))

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',type=open,required=True)
	parser.add_argument('-ar','--ref_annotation',help='annotation  ref file',dest='annoref',type=open,required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	parser.add_argument('-c','--col',help='which col is key, first is input file , second is ref annotation',dest='col',type=int,nargs='+')
	args=parser.parse_args()

	annotation = {}
	annotation , header , ar_len = read_annotation(args.annoref , args.col[1] , 'ref' , annotation)
	for count , line in enumerate(args.input):
		if count == 0 : 
			args.output.write('{0}\t{1}\n'.format(line.rstrip(),header))
		else : 
			tmp = line.rstrip().split('\t')
			name = tmp[args.col[0]]
			# anno ='--\t'*ar_len
			# anno = anno.rstrip()
			if name in annotation:
				anno = annotation[name]
				args.output.write('{0}\t{1}\n'.format(line.rstrip(),anno))
			else:
				print("{0} not in your -ar file, please check".format(name))

if __name__=='__main__':
	main()
