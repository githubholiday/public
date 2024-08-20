#-*- encoding:utf-8 -*-
'''

'''
import argparse
import sys
import os
import time
from datetime import datetime
import glob
bindir = os.path.abspath(os.path.dirname(__file__))

__author__='chengfangtu'
__mail__= 'chengfangtu@genome.cn'

def get_relation_dict( relation_file ):
    relation_dict = {}
    with open( relation_file, 'r') as infile:
        for line in infile:
            tmp = line.rstrip().split('\t')
            map_id = tmp[0].replace("map","")
            if map_id not in relation_dict:
                relation_dict[map_id] = '\t'.join(tmp[1:4])
            else:
                print("{0} is repeat in relation_dict".format(map_id))
    return relation_dict 
def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument('-i','--infile',help='infile of kegg result',dest='infile')
    parser.add_argument('-r','--relation',help='relation of kegg pathway',dest='relation')
    parser.add_argument('-o','--outfile',help='outfile of kegg pathway annotation',dest='outfile')
    args=parser.parse_args()
    relation_dict = get_relation_dict( args.relation )
    with open(args.infile,'r') as input, open(args.outfile, 'w') as output:
        for line in input:
            tmp = line.rstrip().split('\t')
            del tmp[-1]
            if line.startswith("Gene_ID"):
                head = tmp+["Level3","Level2","Level1"]
                output.write('\t'.join(head)+'\n')
                continue
            map_id = tmp[2]
            if map_id in relation_dict:
                level_info = relation_dict
                content = '\t'.join(tmp)+level_info
                output.write(content+'\n')

if __name__ == '__main__':
    main()


