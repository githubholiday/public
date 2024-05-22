#/annogene/datayw/share/software/install/Python-3.3.2/bin/python3
#二代芯片拆分时长统计&二代项目过滤时长统计
'''
v2:增加二代二次拆分的删除和让步监控    zxh at 20230705
'''

import pymysql # 链接mysql，并获取其中的数据库
import argparse # 程序定义它需要的参数，从sys.argv解析出参数
import sys # 获取命令行参数
import os            #应用到系统
import os
import time
import sqlite3
import re
import subprocess
import datetime
import logging
import argparse
import sys

bindir = os.path.dirname(os.path.abspath( __file__ ))
sys.path.append('{0}/../lib'.format(bindir))
#import robot

class Lims():
    def __init__(self):
        #my_limgs = LIMS(config_file)
        usr = 'raw'
        pwd = 'K8b#xx2^QcYs'
        port = 3306
        #host = 'rm-2zeb5kq11g3u77b5bo.mysql.rds.aliyuncs.com'
        host = 'rm-2zeb5kq11g3u77b5b.mysql.rds.aliyuncs.com'
        dbname = 'lims3'
        charset = 'utf8'
        self.db = pymysql.connect(user=usr, password=pwd, host=host, port=port, database=dbname,
                                          charset=charset)
        self.cursor = self.db.cursor()
        self.table1 = 'tb_arrange_toquality'                                                         #集群拆分重拆记录
        self.table2 = 'tb_analysis_task'                            #PM重拆记录
        self.table3 = 'tb_filter_task'                              #过滤记录表

    def get_split_fc(self):
        '''
        选择二代拆分任务
        '''
        #选择二代拆分任务
        cmd1 = "select fc_no,split_start_time,companyId,id,data_dir,generations from {0} where split_end_time IS NULL and split_start_time IS NOT NULL and generations = 0".format(self.table1)
        self.cursor.execute(cmd1)
        data1 = self.cursor.fetchall()

        #选择三代拆分任务
        cmd2 = "select fc_no,split_start_time,companyId,id,data_dir,generations from {0} where split_end_time IS NULL and split_start_time IS NOT NULL and generations = 1".format(self.table1)
        self.cursor.execute(cmd2)
        data2 = self.cursor.fetchall()
        data = data1 + data2
        return data

    def get_resplit_fc( self ):
        cmd3 = "select fc_no,start_time,location,id,type from {0} where start_time IS NOT NULL and type IN ('split','concession','delete') and state IN (2,3)".format(self.table2)     #2是任务投递上了也就是work.sh运行完成，3是任务投递失败也就是work.sh运行失败
        self.cursor.execute(cmd3)
        data3 = self.cursor.fetchall()
        return data3

    def get_filter_pro(self):
        cmd1 = "select merge_project_code, create_time, place, fc_no, id, analysis_path, generations  from {0} where filter_status != 'FILTER_FINISHED' and generations = '2th' and id > 1770".format(self.table3)          #用create_time，即插入起开始计算时间，有些在排队的没有start_time
        self.cursor.execute(cmd1)
        data1 = self.cursor.fetchall()
        cmd2 = "select merge_project_code, create_time, place, fc_no, id, analysis_path,generations  from {0} where filter_status != 'FILTER_FINISHED' and generations = '3th' and id > 1770".format(self.table3)          #用create_time，即插入起开始计算时间，有些在排队的没有start_time
        self.cursor.execute(cmd2)
        data2 = self.cursor.fetchall()

        return data1+data2

    def get_un_filter_pro(self, unfinish_pro):                   #有些id没有核销
        cmd = "select project_code, fc_no from {0} where filter_status = 'FILTER_FINISHED' and generations = '2th' and id > 1770".format(self.table3)              #id1770之前的是旧流程的项目，没有时间和状态的记录
        self.cursor.execute(cmd)
        data_all = self.cursor.fetchall()
        data_un = []
        for i in unfinish_pro:
            #如果未完成的任务不在
            if (i[0], i[3]) not in data_all:
                data_un.append(i)
        return data_un


    def close(self):
        self.cursor.close()
        self.db.close()

class Time():
    def __init__(self):
        self.today = datetime.datetime.now()  # 显示当前时间
        self.now_date = self.today.strftime("%Y-%m-%d %H:%M:%S")  # 将当前时间按一定格式显示成字符串
        self.now_date_today = self.today.strftime("%Y%m%d")

    def difftime(self,starttime):
        diff = str(self.today - starttime)
        if len(diff.split(',')) == 2:
            d = diff.split(',')[0].split()[0]
            h = diff.split(',')[1].split(':')[0]
        else:
            d = 0
            h = diff.split(':')[0]
        return diff, d, h

def print_content( content_list, head):
    print("\t".join(head))
    for tmp in content_list:
        tt = [str(i) for i in tmp ]
        print('\t'.join(tt))

def main():
    aa = Lims()
    #获取所有拆分未完成任务
    split_task = aa.get_split_fc()
    resplit_task = aa.get_resplit_fc()
    #获取所有过滤未完成任务
    filter_task = aa.get_filter_pro()
    u = aa.get_un_filter_pro(filter_task)

    t = Time()
    content = []
    label_list = ['type', 'fc/pro', 'start_time', 'Time_Interval', 'location', 'ID']
    title = 'NGS_split时长统计'
    for i in split_task:
        ttype = "拆分"
        fc = i[0]
        start_time = i[1].strftime("%Y-%m-%d %H:%M:%S")
        #start_time = i[1]
        location = i[2]
        id = i[3]
        path = i[4]
        generation = i[5]
        if generation == 0:
            s_type = "二代"
        else:
            s_type = "三代"
        diff, day, hour = t.difftime(i[1])
        if int(day) < 7:
            #print('split', i[0], i[1], diff.split('.')[0], i[2], i[3])
            content.append([ttype, fc, start_time, diff.split('.')[0], location, id])#, path])
    
    #print_content(content,head)
    resplit_content = []
    for i in resplit_task:
        ttype = "重拆"
        fc = i[0]
        start_time = i[1].strftime("%Y-%m-%d %H:%M:%S")
        location = i[2]
        id = i[3]
        #path = i[4]
        diff, day, hour = t.difftime(i[1])
        if int(day) < 7:
            #print('split', i[0], i[1], diff.split('.')[0], i[2], i[3])
            content.append([ttype, fc, start_time, diff.split('.')[0], location, id])
    print("#拆分和过滤时长超过 【4】小时的需要关注")
    head = ["拆分类型","芯片号","开始拆分时间","拆分时长","地点","ID"]        
    print_content(content,head)
    print("\n")
    #for i in u:
    filter_content = []
    for i in filter_task:
        create_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(i[1]/1000))
        create_time_date = datetime.datetime.strptime(create_time, "%Y-%m-%d %H:%M:%S")
        diff, day, hour = t.difftime(create_time_date)
        generation = i[6]
        f_type= "2th"
        if generation == '3th':
            f_type = "3th"
        if int(day) < 15:
            #print('filter', i[0], create_time_date, diff.split('.')[0], i[2], i[4])
            filter_content.append([f_type, i[0], create_time_date, diff.split('.')[0],i[2],i[4]])#,i[5]])
    head = ["过滤类型","项目号","开始过滤时间","过滤时长","地点","ID"]
    print_content(filter_content, head)
    #robot.Robots(title, label_list, content)
if __name__ == '__main__':
    main()
