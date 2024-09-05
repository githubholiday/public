#-*- encoding:utf-8 -*-
'''
notice文件中输出的 项目过滤目录下没有符合条件的删除记录:均不符合目录的开始过滤时间早于启动交付时间，一般这种是本地目录已经删除。
所以根据给定的项目列表，将tb_delte_info表中的这些项目的local_delete_bool修正为1
'''

import argparse
import sys
import os
import time
from datetime import datetime
import glob
bindir = os.path.abspath(os.path.dirname(__file__))
sys.path.append( '{0}/../../lib'.format( bindir ))
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
    parser.add_argument('-c','--pipe_config',help='config_file',dest='pipe_config',default = '{0}/../../config/config.txt'.format(bindir))
    parser.add_argument('-i','--infile',help='infile',dest='infile')
    args=parser.parse_args()

    my_log.info("配置文件:"+args.pipe_config)
    #读取config文件以及数据库连接
    my_log.info("连接lims数据库")
    lims_db = DOMysql.SQL( args.pipe_config, db_type='lims' )
    my_log.info("连接数据删除数据库")
    local_sql = DOMysql.SQL( args.pipe_config )

    

    #1. 获取tb_deletion_info中的项目信息,交付方式为others的不要
    tb_name = lims_db.config_dic["table"]

    project_list_file = open(args.infile, 'r')
    with open(args.infile, 'r') as input:
        for line in input:
            tmp = line.rstrip()
            update_condition = [('project_id',tmp)]
            update_content = [('local_delete_bool',1)]
            local_sql.update(tb_name, update_content, update_condition)

if __name__ == '__main__':
    main()
