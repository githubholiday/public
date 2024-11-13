'''
根据指定的文件 运行anno go 和 kegg分析

config文件中必须有：
REF_ANNOTATION=注释文件
KEGG_annotate=ko.list
GO_annotate=go.list
relation=gene_id和gene_name对应关系表
category=animal plant fungi
'''

#! /usr/bin/env python3
import argparse
import time
import sys
import re
import os
bindir = os.path.abspath(os.path.dirname(__file__))

__author__ = 'Holiday'
__mail__ = 'Holiday@genome.cn'
__doc__ = 'the description of program'
'''
The up and down gene count in second_level_go 
'''
pat1=re.compile('^s+$')

def info_from_infile( infile):
	'''
	infile:string,输入文件路径
	'''
	infile_path = os.path.abspath(infile)
	infile_dir_path = os.path.dirname(infile_path)
	infile_name = os.path.basename(infile_path).rstrip(".xls")
	#cmp_name = infile_name.split("_gene_symbol.anno.xls")[0]
	#cmp_name = infile_name.split('_diff_gene_symbol.xls')[0]
	return infile_dir_path, infile_name
class Config():
	def __init__(self, config_file ):
		self.config_file = config_file
	def read_config( self ):
		'''
		读取物种的注释配置文件
		获取 ko.list go.list anno_file relation_file category
		'''
		config_dict = {}
		with open( self.config_file, 'r') as infile:
			for line in infile:
				if line == "\n": continue
				tmp = line.rstrip().split('=',1)
				key = tmp[0]
				value = tmp[1]
				if key == "REF_ANNOTATION":
					self.anno = value
				if key == "KEGG_annotate":
					
					self.ko = value
				if key == "GO_annotate":
					self.go = value
				if key == "relation":
					self.relation = value
				if key == "category":
					self.category = value
							
	
def para_check(args, map_dic):
	check_status = 0 

	if len(args.input) == 0 :
		print("输入文件不存在，请检查，退出~")
		check_status += 1
	
	if args.action not in map_dic :
		print("-a must be go or kegg，please retry")
		check_status += 1
	if check_status > 0 :
		print("check your parameters and try again,bye bye~")
		sys.exit(1)

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-i','--input',help='input file',dest='input',nargs='+',required=True)
	parser.add_argument('-o','--outdir',help='outdir of the output',dest='output',required=True)
	parser.add_argument('-s','--shell',help='shell dir',dest='shell',required=True)
	parser.add_argument('--all','--all',help='是否运行',dest='all',action="store_true")
	parser.add_argument('-c','--config',help='category name',dest='config',required=True)
	args=parser.parse_args()
	#是否运行Up和down
	if args.all : 
		type_list = ['all', "Up", "Down"]
	else:
		type_list = ['all']
	#获取物种配置文件中中的信息
	my_conf = Config( args.config )	
	#定义shell文件

	fun_shell = "{0}/function.sh".format(args.shell)
	anno_shell = "{0}/anno.sh".format(args.shell)
	
	fun_shell_handle = open(fun_shell, 'w')
	anno_shell_handle = open(anno_shell, 'w')
	
	for infile in args.input:
		infile_dir_path, file_name = info_from_infile( infile ) #输入文件路径，输入文件名称，
		
		anno_cmd = "make -f {BIN}/anno.mk log_file=LOGFILE sample={sample} infile={infile} outdir={outdir}/{sample}/ relation_file={config.relation} anno_file={config.anno} Anno && echo {sample} Anno finished ".format(BIN=bindir, sample=file_name, infile=infile, outdir=args.outdir,  config=my_conf )
		anno_shell.handle.write( anno_cmd+'\n')

		for up_down in type_list:
			if up_down == "all":
				function_cmd ='''
make -f {BIN}/anno.mk log_file=LOGFILE gene_list={outdir}/{sample}/{sample}.xls go_dir={outdir}/{sample}/{type}/GO sample={sample} go={config.go} relation_file={config.relation} GO && echo {sample} GO finished
make -f {BIN}/anno.mk log_file=LOGFILE gene_list={outdir}/{sample}/{sample}.xls kegg_dir={outdir}/{sample}/{type}/GO sample={sample} ko={config.ko} relation_file={config.relation} category={config.category} KEGG && echo {sample} KEGG finished'''.format(BIN=bindir, sample=file_name, infile=infile, outdir=args.outdir, type=up_down, config=my_conf ) 
			else:
				function_cmd ='''
make -f {BIN}/anno.mk log_file=LOGFILE outdir={outdir}/{sample}/{type} up_or_down={type} infile={outdir}/{sample}/{sample}.xls GetList && echo {sample} gene_list finished
make -f {BIN}/anno.mk log_file=LOGFILE gene_list={outdir}/{sample}/{type}/{type}.gene.xls go_dir={outdir}/{sample}/{type}/GO sample={sample} go={config.go} relation_file={config.relation} GO && echo {sample} GO finished
make -f {BIN}/anno.mk log_file=LOGFILE gene_list={outdir}/{sample}/{type}/{type}.gene.xls kegg_dir={outdir}/{sample}/{type}/GO sample={sample} ko={config.ko} relation_file={config.relation} category={config.category} KEGG && echo {sample} KEGG finished'''.format(BIN=bindir, sample=file_name, infile=infile, outdir=args.outdir, type=up_down, config=my_conf )
			fun_shell_handle.write(function_cmd )
		fun_shell_handle.close()
		anno_shell_handle.close()
		print("输出文件分别为:\n {0}\n{1}".format( anno_shell, fun_shell))
if __name__ == '__main__':
	main()

