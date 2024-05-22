'''
功能：
1. 根据给定的芯片号统计所有芯片的产出情况，以及reads的产出得率
2. 统计tb_arrange_lane数据库中所有二代芯片的产出情况
需要的数据库:
tb_arrange_lane 


'''
import os
import sys
import argparse
import glob
import time
bindir = os.path.dirname(os.path.abspath( __file__ ))
filename=os.path.basename(__file__)
from Lims_SQL import LIMS

file_name = os.path.basename(__file__)

__author__ = 'Tu chengfang'
__mail__ = 'chengfangtu@genome.cn'

def get_fc_list( fc_list_file ):
    fc_list = []
    with open(fc_list_file, 'r') as infile:
        for line in infile:
            tmp = line.rstrip()
            fc_list.append(tmp)
    return fc_list

def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument( '-f', '--flist', help='fc list file, like HCSF*', dest = 'flist')
    parser.add_argument( '-c','--config',help='config',dest='config',default="{0}/config.ini".format(bindir) )
    parser.add_argument( '-o','--output',help='output',dest = 'output',required=True )
    args = parser.parse_args()


    lims_db = LIMS( config_file = args.config )
    
    select_list = ["fc_no","lane","task_state","lib_num","undetermined_reads","lib_product_reads","total_reads","phix_reads","undetermined_base","lib_product_base","total_base","phix_base"]
    #从中间表中获取芯片信息，决定是否运行后面的监控步骤
    if args.flist:
        all_records = []
        fc_list = get_fc_list( args.flist )
        for fc in fc_list:
            conditions = [("fc_no",fc)]
            records = lims_db.select(table_name="tb_arrange_lane", col_list=select_list, conditions=conditions)
            if len(records) == 0 :
                print("{0} not in tb_arrange_lane".format( fc ))
            else:
                [ all_records.append(i) for i in records ]
    else:
        all_records = lims_db.select( table_name='tb_arrange_lane',col_list=select_list, conditions = None)
    ["fc_no","lane","task_state","lib_num","undetermined_reads","lib_product_reads","total_reads","phix_reads","undetermined_base","lib_product_base","total_base","phix_base"]
    with open( args.output, 'w') as outfile:
        head = ["fc_no","lane","task_state","lib_num","undetermined_reads","lib_product_reads","total_reads","phix_reads","undetermined_base","lib_product_base","total_base","phix_base","Rate"]
        outfile.write("\t".join(head)+'\n')
        for tmp in all_records:
            if not tmp[5] : continue #三代芯片没有这几个数据
            lib_product_reads = float(tmp[5])
            total_reads = float(tmp[6])
            tt = [str(i) for i in tmp ]
            rate = '{0:.2f}'.format((lib_product_reads/total_reads)*100)
            outfile.write("\t".join(tt)+'\t'+rate+'\n')


                                                                                                                                                   
        
if __name__ == '__main__':
    main()
