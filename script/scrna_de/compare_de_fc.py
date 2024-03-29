
'''
功能：
'''


import argparse
import sys # 获取命令行参数
import os
import re



bindir = os.path.dirname(os.path.abspath( __file__ ))
sys.path.append('{0}/../lib'.format(bindir))
#import robot
def get_gene_dic( gene_file ):
    gene_dict = {}
    with open( gene_file, 'r') as infile:
        for line in infile:
            if line.startswith(','): continue
            tmp = line.rstrip().split(',')
            gene_name = tmp[0]
            fc = float(tmp[1])
            abs_fc = abs(fc)
            if gene_name not in gene_dict:
                gene_dict[gene_name] = abs_fc
    return gene_dict

def compare_gene( dict1, dict2, recue_file, other_file ):
    with open( recue_file, 'w') as out1, open( other_file, 'w') as out2:
        out1.write('GeneName')
        for gene_name in dict1:
            if gene_name not in dict2:
                out2.write(gene_name+'\n')
                continue
            lg_fc1 = dict1[gene_name]
            lg_fc2 = dict2[gene_name]
            if lg_fc1 > lg_fc2 :
                out1.write(gene_name+'\n')

def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument( '-i1', '--infile1', help='de gene file of cmp1', dest='infile1', required=True)
    parser.add_argument( '-i2', '--infile2', help='de gene file of cmp2', dest='infile2', required=True)
    parser.add_argument( '-o', '--outdir', help='outdir of result', dest='outdir', required=True)
    parser.add_argument( '-p', '--prefix', help='outdir of result', dest='prefix', required=True)
    args = parser.parse_args()

    de_gene_dict1 = get_gene_dic( args.infile1 )
    de_gene_dict2 = get_gene_dic( args.infile2 )
    rescue_file = '{0}/{1}_rescue.txt'.format( args.outdir, args.prefix)
    #在文件1中，但是不在文件2中
    not_in_dict2_file = '{0}/{1}_not_in_2.txt'.format( args.outdir, args.prefix)
    compare_gene(de_gene_dict1,de_gene_dict2, rescue_file,not_in_dict2_file)





if __name__ == '__main__':
    main()
