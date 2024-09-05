#-*- encoding:utf-8 -*-
'''
由于过滤流程bug，导致样本数量合格就上云（即样本量合格，不考虑数据量是否合格）

该脚本是从删除数据库中获取项目列表，并从Filter_Result/email/email_feedback.txt中获取finish是否为yes
如果是yes，则跳过
如果是no,则输出列表

下一步应该是删除这些输出列表中的信息
'''

import argparse
import sys
import os
import time
from datetime import datetime
import glob
bindir = os.path.abspath(os.path.dirname(__file__))
sys.path.append( '{0}/../lib'.format( bindir ))
import DOMysql
from logger import Log
import robot

my_log = Log( 'delete_fq.log' )

__author__='chengfangtu'
__mail__= 'chengfangtu@genome.cn'

    

def read_email_feedback( infile ):
    yes_or_no = ''
    with open( infile, 'r') as input:
        for line in input:
            if line.startswith("finish"):
                tmp = line.rstrip().split('=')
                yes_or_no = tmp[1]
    return yes_or_no
def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument('-c','--pipe_config',help='config_file',dest='pipe_config',default = '{0}/../config/config.txt'.format(bindir))
    parser.add_argument('-o','--outfile',help='outfile',dest='outfile')
    args=parser.parse_args()

    my_log.info("配置文件:"+args.pipe_config)
    #读取config文件以及数据库连接
    my_log.info("连接lims数据库")
    lims_db = DOMysql.SQL( args.pipe_config, db_type='lims' )
    my_log.info("连接数据删除数据库")
    local_sql = DOMysql.SQL( args.pipe_config )

    output = open(args.outfile, 'w')

    #1. 获取tb_deletion_info中的项目信息,交付方式为others的不要
    tb_name = lims_db.config_dic["table"]
    select_condition1 = [("delivery_account","ali")]
    select_condition2 = [("delivery_account","hw")]
    select_condition = select_condition1+select_condition2
    records1 = local_sql.select( tb_name, "*",select_condition1 )
    records2 = local_sql.select( tb_name, "*",select_condition2 )
    records = records1 + records2 


    for tmp in records :
        project_id = tmp[1]
        cloud_address = tmp[3]
        delivery_dir_name = cloud_address.rstrip('/').split('/')[-1]
        delivery_dir_name_split = delivery_dir_name.split("_")
        #2. 过滤数据：2）cloud_address中的最后目录使用_分割，length超过3的不要；
        if len(delivery_dir_name_split) > 3 : continue
        #3. 根据筛选的项目从tb_filter_task中获取项目过滤目录，并获取其Filter_Result目录
        filter_select_condition = [("project_code", project_id)]
        filter_records = lims_db.select("tb_filter_task", ["analysis_path"], filter_select_condition)
        if len(filter_records) == 0 :
            print("tb_filter_task 没有 {0} 信息".format( project_id ))
            continue
        analysis_path = filter_records[0][0]
        filter_result_path = "{0}/../Filter_Result".format(analysis_path)
        if not os.path.exists(filter_result_path) :
            print("{0} 路径不存在".format( filter_result_path ))
            continue
        email_feedbakc_file = '{0}/email/email_feedback.txt'.format(filter_result_path )
        if not os.path.exists(email_feedbakc_file):
            print("{0} 路径不存在".format( email_feedbakc_file ))
            continue
    #4. 读取Filter_Result/email/email_feedback.txt中DATA,finish=
        finish_status = read_email_feedback(email_feedbakc_file)
        if finish_status == "YES": continue
    #finish=yes，如果是yes，则可以删除，如果不是，输入id，项目编号，等信息
        del_condition = [("ID",tmp[0])]
        print(del_condition)
        local_sql.delete( tb_name, del_condition)
        out_value = [str(i) for i in tmp]
        output.write("\t".join(out_value)+'\n')
    #4. 读取Filter_Result/email/email_feedback.txt中DATA,finish=yes，如果是yes，则可以删除，如果不是，输入id，项目编号，等信息，以及Filter_Result/email/email_feedback.txt的修改时间
    #最后出来的表格需要核对，核对之后才能确认哪些id是可以删除的

    output.close()

if __name__ == '__main__':
    main()



