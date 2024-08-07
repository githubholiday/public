'''
输入samtools的flagstat文件和 bedtools的bamcoverage文件，进行汇总
flagstat:
9721892 + 0 in total (QC-passed reads + QC-failed reads)
9721892 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
111973 + 0 mapped (1.15% : N/A)
111973 + 0 primary mapped (1.15% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

depth:
genome  0       14476   139962  0.103428
genome  1       407     139962  0.00290793
genome  2       355     139962  0.0025364
genome  3       267     139962  0.00190766
genome  4       436     139962  0.00311513


'''


#! /usr/bin/env python3
import argparse
import sys
import os
import re
bindir = os.path.abspath(os.path.dirname(__file__))

__author__='Liu Tao'
__mail__= 'liutao@seqwisdom.com'

pat1=re.compile('^\s*$')


def parse_flagstat(flagstat):
	total_reads = 0
	mapped_reads = 0
	for line in flagstat:
		if pat1.match(line):
			continue
		if 'in total' in line:
			total_reads = int(line.split()[0])
		if 'primary mapped' in line:
			mapped_reads = int(line.split()[0])
	return total_reads, mapped_reads

def parse_depth(infile , lower_than ):
	total_genome = 0
	lower_than_count = 0
	depth = open(infile)
	for line in depth:
		if pat1.match(line):
			continue
		if line.startswith('genome'):
			tmp = line.split()
			total_genome = int(tmp[3])
			base_depth = int(tmp[1])
			if base_depth <= lower_than:
				lower_than_count += int(tmp[2])
	depth.close()
	return total_genome, lower_than_count


def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-f','--flagstat',help='input flagstat',dest='flagstat',type=open,required=True)
	parser.add_argument('-d','--depth', help='input depth',dest='depth',required=True)
	parser.add_argument('-n','--name',help='input name',dest='name',type=str,required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	#parser.add_argument('-m','--mm',help='output file',dest='mm',action='store_false')
	args=parser.parse_args()

	total_reads , mapped_reads = parse_flagstat(args.flagstat)
	mapped_ratio = mapped_reads / total_reads * 100

	total_genome, zero_depth = parse_depth(args.depth , 0 )
	coverage_base = total_genome - zero_depth
	coverage = (total_genome - zero_depth) / total_genome * 100

	total_genome, lower_than_4_depth = parse_depth(args.depth , 4 )
	higher_than_4_coverage_genome  = total_genome - lower_than_4_depth
	higher_than_4_coverage_ratio  = (total_genome - zero_depth) / total_genome * 100

	args.output.write('Name\t{0}\n'.format(args.name))
	args.output.write('Total_Reads\t{0}\n'.format(total_reads))
	args.output.write('Mapped_Reads\t{0}\n'.format(mapped_reads))
	args.output.write('Mapped_Ratio\t{0:.2f}\n'.format(mapped_ratio))
	args.output.write('Total_Genome\t{0}\n'.format(total_genome))
	args.output.write('Coverage_Genome\t{0}\n'.format(coverage_base))
	args.output.write('Coverage_Ratio\t{0:.2f}\n'.format(coverage))
	args.output.write('Coverage_Genome(>4)\t{0}\n'.format(higher_than_4_coverage_genome))
	args.output.write('Coverage_Ratio(>4)\t{0:.2f}\n'.format(higher_than_4_coverage_ratio))





	

if __name__ == '__main__':
	main()
