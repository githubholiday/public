'''
获取素有基因的kegg信息，输出两个文件，分别
1）基因和其对应的K信息
2）基因和其对应的通路，map信息
'''

import argparse

__author__ = "chengfangtu"
__mail__ = "ruiliao@genome.cn"

def deal_kegg( infile, ko_handle, pathway_handle):
	pathway_value = []
	flag = 0
	gene_name = ''
	K_id = ''
	with open(infile, 'r') as infile :
		for line in infile :
			if line.startswith('#') : continue
			if line.startswith('gene') :
				tmp = line.rstrip().split('\t')
				gene_name = tmp[0]
				k_info = tmp[1]
				if k_info == 'None' : continue
				if not k_info.startswith('K') : continue
				K_id = k_info.split('|')[0]
				ko_value = [ gene_name, K_id]
				ko_handle.write('\t'.join(ko_value)+'\n')
			elif line.startswith('///'):
				flag = 0
				gene_name = ''
				K_id = ''
			elif line.startswith('Query'):
				gene_name = line.rstrip().split('\t')[1]
				pathway_value = []
				pathway_value.append(gene_name)
				flag = 1
			elif line.startswith('KO:') and flag == 1:
				K_id = line.rstrip().split('\t')[1]
			elif line.startswith('Pathway') and flag == 1 :
				tmp = line.rstrip().split('\t')
				kegg_description = tmp[1]
				map_id = tmp[3].replace('ko','')
				pathway_value = [ gene_name, K_id, map_id,  kegg_description]
				pathway_handle.write('\t'.join(pathway_value)+'\n')
			elif flag == 1 :
				tmp = line.rstrip().split('\t')
				kegg_description = tmp[1]
				map_id = tmp[3].replace('ko','')
				pathway_value = [ gene_name, K_id, map_id,  kegg_description]
				pathway_handle.write('\t'.join(pathway_value)+'\n')
				


def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='the input file of kobas out',dest='input',nargs='+',required=True)
	parser.add_argument('-o','--output',help='pathway file',dest='out',type=argparse.FileType('w'),required=False)
	parser.add_argument('-k','--koout',help='kofile file',dest='ko',type=argparse.FileType('w'),required=False)
	args=parser.parse_args()

	path,ko = [],[]
	args.ko.write("Gene_ID\tAnnotation\n")
	args.out.write("Gene_ID\tAnnotation\tMap\tPathway\n")
	for infile in args.input:
		deal_kegg( infile, args.ko, args.out)
'''
		dic1,dic2 = read_file(infile, args.ko, args.out)
		path.extend(dic1)
		ko.extend(dic2)
	if args.ko:
		args.ko.write("Gene_ID\tAnnotation\n")
		ko_al = []
		for i in ko:
			if i not in ko_al:
				args.ko.write("\t".join(i)+"\n")
				ko_al.append(i)
	if args.out:
		args.out.write("Gene_ID\tAnnotation\tMap\tPathway\n")
		path_dic = []
		for i in path:
			if i not in path_dic:
				args.out.write("\t".join(i)+"\n")
				path_dic.append(i)
'''
if __name__ == "__main__":
	main()
