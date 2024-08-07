import argparse

__author__ = "liaorui"
__mail__ = "ruiliao@genome.cn"

def read_file(file1,nu):
	gene = []
	with open(file1) as infile:
		for line in infile:
			tmp = line.strip().split("\t")
			if int(tmp[1]) > nu:
				gene.append(tmp[0])
	return gene

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',nargs='+',required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	parser.add_argument('-n','--num',help='count num',dest='nu',type=int,required=False,default=100)
	args=parser.parse_args()

	gene = []
	for name in args.input:
		gene.extend(read_file(name,args.nu))
	gene_al = list(set(gene))
	for i in gene_al:
		args.output.write(i+"\n")

if __name__ == "__main__":
	main()
