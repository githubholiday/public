#-*- encoding:utf-8 -*-
'''
【功能】
1. 获取删除数据库中交付未完成的项目，并输出启动交付时间超过7天的，输出项目编号，云上交付地址，交付数据库id，启动交付时间，是否交付完成, 交付完成时间
2. 根据上述信息，人工进行check，并将最后一列修改为实际情况，yes or no ,根据该表更新数据删除数据库中的delivery_bool=1，delivery_finish_time=实际时间，按照 2024-5-30 14:40:30格式填写


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

my_log = Log( 'get_exception.py' )

__author__='chengfangtu'
__mail__= 'chengfangtu@genome.cn'

class TimeFormat():
    '''
    功能：时间的模块,可以获取时间戳、目录创建时间等
    参数：
    date_format:时间格式，默认为'%Y-%m-%d'
    '''
    def __init__(self, date_format = '%Y-%m-%d %H:%M:%S'):
        self.date_format = date_format
    def now_time(self):
        '''
        功能：获取当前时间
        return:当前时间，str类型,格式为2024-05-13
        '''
        return datetime.now().strftime(self.date_format)
    def get_delta_days( self, time1, time2 ):
        '''
        功能：计算两个时间(时间2-时间1)之间的天数差值
        参数：
            time1:时间1，str类型，格式为2024-05-13-09-09-09
            time2:时间2，str类型，格式为2024-05-13-
        return:时间差值，int类型
        '''
        delta = datetime.strptime(time2, self.date_format) - datetime.strptime(time1, self.date_format)
        return delta.days
    def get_delta_seconds( self, time1, time2 ):
        delta = datetime.strptime(time2, self.date_format) - datetime.strptime(time1, self.date_format)
        delta_seconds = delta.seconds + delta.days*24*3600
        return delta_seconds
    def dir_create_time( self, dir_path ):
        '''
        功能：获取一个目录的创建2024-05-13
        参数：
            dir_path:目录路径，str类型，格式为/home/chengfangtu/test
        返回值：
            创建时间，str类型，格式为2024-05-13
        '''
        create_stamptime = os.path.getmtime(dir_path) #时间戳 163467344
        create_localtime= time.localtime(create_stamptime/1000)
        create_formatime = time.strftime(self.date_format ,create_localtime) #格式化为时间 2024-04-14
        return create_formatime

class SelectItem():
    def __init__(self, outfile, local_sql, config_dict, date_format = '%Y-%m-%d %Y-%m-%d'):
        self.local_sql = local_sql
        self.date_format = date_format
        self.config_dict = config_dict
        self.outfile = outfile
        self.time_handle = TimeFormat(date_format)
    def main( self ):
        '''
        功能：获取未交付完成，且启动交付时间超过7天的项目，并输出到self.outfile中
        返回值：None
        '''
        all_records = self.get_all_records()
        flag = 0 
        with open(self.outfile, 'w') as output:
            head = ["#删除数据库ID","子项目编号","云上交付地址","交付数据库ID","启动交付时间","是否交付完成","交付完成时间"]
            output.write('\t'.join(head)+'\n')
            for tmp in all_records:
                id = tmp[0]
                project_id = tmp[1]
                delivery_start_time = tmp[2] #启动交付时间
                cloud_address = tmp[3]
                delivery_db_id = tmp[7]
                delivery_end_time = tmp[9] if tmp[9] else "" 
                now_time = self.time_handle.now_time()
            #获取当前时间-启动交付时间的delta天数
                self.delta_days = self.time_handle.get_delta_days( delivery_start_time,now_time)# now_time)
                #如果是超期删除，则将该项目的analysis_bool是否为分析项目的值改成0，即非分析项目，后面按照非分析项目的删除规则执行判断
                if self.delta_days > 7:
                    out_value = [str(id), project_id, cloud_address, str(delivery_db_id),delivery_start_time, "no", "" ]
                    output.write('\t'.join(out_value) + '\n')
                    flag =1 
            if flag == 0 :
                my_log.info("没有交付7天仍未完成的项目")

    def get_all_records(self):
        '''
        功能：获取数据库中 交付是否完成=否
        return:所有项目列表，list类型[(),(),()]
        '''
        select_conditions = [ ("delivery_bool",0) ]
        tb_name = self.config_dict["table"]
        all_records = self.local_sql.select( table_name=tb_name, conditions=select_conditions )
        if len(all_records ) == 0 :
            my_log.info( "{0}中没有交付未完成的条目" )
            return []
        else:
            return all_records

                      
class Update():
    def __init__( self, infile, local_sql, config_dict, date_format ):
        self.infile = infile
        self.config_dict = config_dict
        self.local_sql = local_sql
        self.time_format = date_format
        self.tb_name = local_sql.config_dic['table']
        self.time_handle = TimeFormat(self.time_format)

    def main(self):
        '''
        功能:根据提取到待删除项目信息（project_dict)，获取删除列表，并进行删除
        过程：
            第一步：根据project_dict中的子项目编号，调用get_filter_dir从lims表tb_filter_task中获取到项目的过滤路径（为Analysis_1/../）
            第二步：获取到过滤路径下的所有文件和目录，调用filter_glob将所有文件和Filter_Result目录过滤掉，只保留目录
            第三步：调用delete_fq处理第二步获取的目录，获取删除列表并将其删除
        参数：
            self.project_dict：待处理项目字典，project_dict={project_id:{"delivery_start_time":"2024-5-15","ID":[1,2,3,4]}
        '''
        with open( self.infile, 'r') as input:
            for line in input:
                tmp = line.rstrip().split('\t')
                if line.startswith("#") : continue
                id = tmp[0]
                delivery_bool = tmp[5]
                delivery_finish_time = tmp[6]
                if delivery_finish_time == '' :
                    delivery_finish_time = self.time_handle.now_time()
                if delivery_bool == 'yes':
                    update_info = [("delivery_bool",1), ("delivery_finish_time", delivery_finish_time)]
                    update_condition = [("ID",id)]
                    self.local_sql.update(table_name= self.tb_name, value_list=update_info, conditions=update_condition)
                #如果确实是交付未完成，则标记为no，删除数据库不更新
                else:
                    continue

def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument('-c','--pipe_config',help='config_file',dest='pipe_config',default = '{0}/../config/config.txt'.format(bindir))
    parser.add_argument('--get',help='get release exception records',dest='get',action='store_true')
    parser.add_argument('--update',help='get release exception records',dest='update',action='store_true')
    parser.add_argument('-i','--infile',help='infile of update',dest='infile')
    args=parser.parse_args()
    if not args.get and not args.update:
        my_log.error("至少提供 --get --update 中一项")
        sys.exit(1)
    
    if args.update and not args.infile :
        my_log.error("需要提供待更新的项目列表文件,使用-i参数")
        sys.exit()
    
    my_log.info("配置文件:"+args.pipe_config)
    #读取config文件以及数据库连接
    my_log.info("连接数据删除数据库")
    local_sql = DOMysql.SQL( args.pipe_config )
    config_dict = local_sql.config_dic
    date_format = config_dict["TimeFormat"]
    #

    if args.get :
        time_handle = TimeFormat( "%Y-%m-%d" )
        now_time = time_handle.now_time()
        outdir = local_sql.config_dic['exception_dir']

        release_exception_file = '{0}/{1}_release_exception.xls'.format( outdir, now_time )

        my_get = SelectItem( release_exception_file, local_sql, config_dict, date_format=date_format)
        my_get.main()
        my_log.info("交付异常项目列表文件: {0}".format(release_exception_file))
        
        robot_url = local_sql.config_dic['robot_url']
        title = "{0} 数据交付异常检测脚本已经执行完毕".format( now_time )
        label = ['日期','脚本',"文件"]
        content = [now_time,os.path.dirname(__file__),release_exception_file ]
        robot.Robots(title=title,label_list=label,content_list=content, url=robot_url)

    if args.update:
        if args.infile :
            my_update = Update( args.infile, local_sql, config_dict, date_format=date_format)
            my_log.info("更新完成")
        else:
            my_log.error("需要提供待更新的项目列表文件")
            sys.exit()

if __name__ == '__main__':
    main()

