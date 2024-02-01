#/annogene/datayw/share/software/install/Python-3.3.2/bin/python3
#二代芯片拆分时长统计&二代项目过滤时长统计
'''
1. 监控北京和义乌 二三代拆分和过滤正在进行的任务，并输出任务运行时长

更新：
2024-1-26 涂成芳
1. 增加三代拆分和过滤的运行时长监控
2. 更改屏幕输出方式，将datetime转化为时间

PENDING_FILTER的显示的时间特别长，有问题
后面加上过滤状态，另外以create_time
'''

#import pymysql # 链接mysql，并获取其中的数据库
import argparse # 程序定义它需要的参数，从sys.argv解析出参数
import sys # 获取命令行参数
import os            #应用到系统
import os
import time
import sqlite3
import re
import subprocess
from datetime import datetime
import logging
import argparse
import sys
from LimsSQL import LIMS
import time

bindir = os.path.dirname(os.path.abspath( __file__ ))
sys.path.append('{0}/../lib'.format(bindir))
#import robot

class Job():
    def __init__(self, config_file ):
        self.my_lims = LIMS( config_file )
        self.time_format = "%Y-%m-%d %H:%M:%S" 

    def get_split_task(self, split_type = '2th'):
        '''
        获取正在拆分中的任务
        '''
        generations = '0'
        if split_type == '3th':
            generations = '1'
        col_list = ["fc_no","split_start_time","split_end_time","companyId","id","data_dir","generations"]
        select_condition = [("taskStatus","=","5"),("generations","=",generations)]
        select_content = self.my_lims.select( 'tb_arrange_toquality', col_list, select_condition )
        if len(select_content) == 0 :
            print("没有正在拆分的 {0} 任务".format( split_type ))
        return select_content

    def get_resplit_task(self):
        '''
        获取重拆表里的任务，之前是只获取split,concession,delete，该版本没有限制，所有都筛选出来
        '''
        select_condition1 = [("state","!=", "4"),("id",">","19660")]
        col_list = ["fc_no","start_time","create_time","entity_id","location","id","type","record"]
        select_content1 = self.my_lims.select( 'tb_analysis_task', col_list, select_condition1 )
        #select_condition2 = [("start_time","IS","NOT NULL"),("state","=", "3")]
        #select_content2 = self.my_lims.select( 'tb_analysis_task', col_list, select_condition2 )
        select_content = select_content1 #+ select_content2
        if len( select_content ) == 0 :
            print("没有正在重新处理的任务")
        return select_content
    
    def get_filter_task(self, filter_type='2th'):
        '''
        获取正在过滤或者过滤异常的任务
        '''
        generations = '2th'
        if filter_type == '3th':
            generations = '3th'
        #cmd1 = "select merge_project_code, create_time, place, fc_no, id, analysis_path, generations,filter_status  from {0} where filter_status != 'FILTER_FINISHED' and generations = '2th' and id > 1770".format(self.table3)          #用create_time，即插入起开始计算时间，有些在排队的没有start_time
        select_condition = [("filter_status","!=","FILTER_FINISHED"),("generations","=", generations),("id",">","1770")]
        col_list = ["merge_project_code","create_time","place","fc_no","id","analysis_path","filter_status"]
        select_content = self.my_lims.select( 'tb_filter_task', col_list, select_condition )
        filter_list = []
        if len(select_content) == 0 :
            print("没有正在重新过滤的任务")
        #再将任务在tb_filter_task搜一遍，查看相同的项目编号和芯片号的记录是否有过滤完成，如果有，则去掉这条记录，如果没有，则继续
        else:
            for select_t in select_content:
                merge_project_code = select_t[0]
                fc_no = select_t[3]
                select_condition = [("merge_project_code","=",merge_project_code),("fc_no","=", fc_no),("filter_status","=","FILTER_FINISHED")]
                project_select = self.my_lims.select( 'tb_filter_task', "*", select_condition )
                if len(project_select) == 0 :
                    filter_list.append(select_t)
        return filter_list
    

    def deal_filter_task( self, filter_task, generation ):
       
       #col_list = ["merge_project_code","create_time","place","fc_no","id","analysis_path"]
       head = ["类型","项目编号","芯片号","开始过滤时间","过滤时长","过滤状态"]
       print("\t".join(head))
       for tmp in filter_task:
           start_time = tmp[1]
           id = str(tmp[4])
           timestamp = time.localtime(start_time/1000)
           delta_hours = delta_time_period( start_time )
           start_time_f = time.strftime(self.time_format,timestamp)
           out_value = [generation, tmp[0], tmp[3], start_time_f, '{0:.2f}'.format(delta_hours),tmp[6],id ]
           print("\t".join(out_value))


    def deal_spit_task(self, split_task_list, generation ):
        '''
        二三代重拆任务的共同处理部分
        '''
        head = ["芯片类型","芯片号","开始拆分时间","拆分时长"]
        print("\t".join(head))
        for tmp in split_task_list:
            fc_no = tmp[0]
            split_start_time = tmp[1]
            split_end_time = tmp[2]
            generations = tmp[6]
            if split_end_time : continue
            start_time = split_start_time.strftime(self.time_format)
            delta_hours = delta_time( split_start_time )
            out_value = [ generation+"    ", fc_no, start_time, '{0:.2f}'.format(delta_hours) ]
            print("\t".join(out_value))


    def deal_2th_split_task( self ):
        '''
        处理二代拆分任务
        '''
        print('\n##二代拆分任务监控')
        ngs_split_task = self.get_split_task(split_type='2th')
        self.deal_spit_task(ngs_split_task,generation='2th')

    def deal_3th_split_task( self ):
        '''
        处理三代拆分任务
        '''
        print('\n##三代拆分任务监控')
        tgs_split_task = self.get_split_task(split_type='3th')
        self.deal_spit_task(tgs_split_task,generation='3th')

    def deal_2th_filter_task(self):
        print('\n##二代过滤任务监控')
        filter_task = self.get_filter_task('2th')
        self.deal_filter_task( filter_task, '2th' )

    def deal_3th_filter_task(self):
        print('\n##三代过滤任务监控')
        filter_task = self.get_filter_task('3th')
        self.deal_filter_task( filter_task, '3th' )
    
    def deal_resplit_task( self ):
        print("\n二三代、让步等处理任务监控")
        resplit_task = self.get_resplit_task()
        col_list = ["fc_no","start_time","entity_id","location","id","type","record"]
        head = ["重拆id","芯片号","项目编号","创建时间","时长","重拆类型","信息"]
        print("\t".join(head))
        for tmp in resplit_task:
            fc_no = tmp[0]
            start_time = tmp[1]
            create_time = tmp[2]
            project_id = tmp[3]
            id = tmp[5] #int类型
            resplit_type = tmp[6]
            record = "" if tmp[7] == None else tmp[7]
            delta_hours = delta_time( create_time )
            start_time_f = create_time.strftime(self.time_format)
            out_value = [str(id), fc_no, project_id, start_time_f, '{0:.2f}'.format(delta_hours), str(resplit_type), record ]
            print("\t".join(out_value))


    def main_job(self, select_type = 'all') :
        '''
        主任务，select_type 可选 all, split, filter
        其中：all表示选择拆分，重拆，过滤任务
        split表示选择 二三代拆分任务
        filter表示选择 二三代过滤任务
        '''
        if select_type == 'split':
            print("#二三代拆分任务如下，关注时长超过 4 小时的任务")
            self.deal_2th_split_task()
            self.deal_3th_split_task()
            self.deal_resplit_task()
        elif select_type == "filter":
            self.deal_2th_filter_task()
            self.deal_3th_filter_task()
        elif select_type == "all":
            self.deal_2th_split_task()
            self.deal_3th_split_task()
            self.deal_resplit_task()
            self.deal_2th_filter_task()
            self.deal_3th_filter_task()
        else:
            print("请输出  split, filter, all 中任意一个")
        
def delta_time( given_time ):
    '''
    计算datetime的时间差
    '''
    now_time = datetime.now()
    delta_seconds = (now_time - given_time).seconds
    delta_hours = delta_seconds/3600
    return delta_hours

def delta_time_period( given_time ):
    '''
    计算毫秒时间戳的时间差，返回hours
    '''
    time_delta = int(time.time()*1000) - given_time
    delta_hours = time_delta/3600000
    return delta_hours


def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument( '-t', '--stype', help='type of select', dest = 'stype',default='all')
    parser.add_argument( '-c','--config',help='config',dest='config',default="{0}/config.txt".format(bindir) )
    args = parser.parse_args()

    my_job = Job( args.config )
    my_job.main_job( args.stype )


if __name__ == '__main__':
    main()
