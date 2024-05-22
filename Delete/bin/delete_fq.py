#-*- encoding:utf-8 -*-
'''
【功能】
根据删除数据库中的信息获取符合条件的项目数据（过滤路径（Analysis和CHIP）和拆分路径），并进行删除

1. 筛选符合条件的项目列表
交付是否完成=是
是否删除=否

2. 计算1中所有项目列表的交付完成时间与当前时间的差值delta
如果是超期删除的项目，则按照非分析项目的判断条件执行
1)如果是分析项目，>30天(analysis_delete_delta),存入待删除字典,
2)如果是非分析项目,>7天(配置文件中的filter_delete_delta),存入待删除字典
=>获取待删除字典 good_records_dict[project_id] = {"delivery_start_time":delivery_start_time, "ID":[1,2,3,4]} 
=>并且启动交付时间获取的最近的一次，启动交付时间为最新一次的未删除记录的启动交付时间

3. 
1）通过2中的待删除列表，从tb_filter_task列表中，获取到该项目的过滤路径 （tb_filter_task.filter_path/../）
2）glob过滤路径下的所有目录，包括Analysis_1和Filter_Result，只筛选Analysis*开头的目录，用于下一步分析
3）获取每个目录的创建时间，如果启动交付时间-创建时间>0即创建在前，启动交付在后(函数：get_delete_list），则获取其对应的芯片号，并从tb_filter_task中获取对应的拆分数据路径和Filter_Result/CHIP/CHIP_FC路径，存入到待删除列表中
4)如果待删除列表为空，则不执行删除操作 *** 需要关注为什么会有不删除的操作
4）执行删除：如果test模式，则print出删除的目录；如果是非test模式，直接os.system('rm -rf 目录')--循环三次（参数）进行删除
5）其中任意一个路径删除失败，该项目为删除失败，并输出到outdir/delete_failed_年-月-日.txt文件中，【可以人工核验删除，按照删除成功更新数据库信息--暂时不处理,将ID也输出】====希望加上机器人提醒
5）删除成功后,根据good_records_dict中的id_list更新删除数据库中该项目的删除状态local_delete_bool=1,delete_time=当前时间(年-月-日) --- 如果有多条记录，多条记录均更新

4. 获取数据库中交付异常的项目 === 单独写一个脚本，不在该脚本中增加
筛选 交付是否完成=否，是否删除=否 的项目列表
判断项目的启动交付时间与当前时间差值，如果大于7天，则输出到outdir/deliver_exception_年-月-日.txt文件中

【参数说明】
-c,--pipe_config:配置文件，使用到配置文件中的信息有 数据库信息,删除时间中的filter_delete_time

【注意事项】
1. 删除项目时，如果同一个项目有多条未删除记录，则以最新的启动交付时间为准
2. 同一个项目，任意一个目录删除失败，则为失败，数据库的删除状态不更新
3. 如果真有删除失败的，可能会导致拆分处的数据没有删除，而且发现不了，暂时搁置，但是
4. 过滤项目重过滤的时候，将Analysis目录删除重新创建
5. .e文件中会输出删除失败的命令，以及没有删除列表的项目编号
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
from lib.logger import Log
import robot

my_log = Log( 'delete_fq.log' )

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
    def delta_days( self, time1, time2 ):
        '''
        功能：计算两个时间(时间2-时间1)之间的天数差值
        参数：
            time1:时间1，str类型，格式为2024-05-13-09-09-09
            time2:时间2，str类型，格式为2024-05-13-
        return:时间差值，int类型
        '''
        delta_seconds = (datetime.strptime(time1, self.date_format) - datetime.strptime(time2, self.date_format)).seconds()
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

class SelectItem():
    def __init__(self, local_sql, config_dict, date_format = '%Y-%m-%d %Y-%m-%d'):
        self.local_sql = local_sql
        self.date_format = date_format
        self.config_dict = config_dict
        self.fiter_delete_delta = config_dict['filter_delete_delta']
        self.analysis_delete_delta = config_dict['analysis_delete_delta']
        self.time_handle = TimeFormat(date_format)
    def main( self ):
        '''
        功能：获取待删除记录的串写
        过程：
            第一步：获取满足条件的所有记录(是否交付完成=是，本地是否删除=否)
            第二步：根据是分析项目和非分析项目的不同判定条件，获取到待删除项目记录，并返回。如果同一个项目都多条记录，则会依次判断启动交付时间，取最近的一次启动交付时间
        返回值:
            待删除项目记录：dict类型-双层字典，格式示例如下：
            good_records_dict[project_id] = {"delivery_start_time":delivery_start_time, "ID":[1,2,3,4]}
        '''
        all_records = self.get_all_records()
        good_records_dict = {}
        for tmp in all_records:
            project_id, delivery_start_time, id  = self.get_good_record(self.fiter_delete_delta, self.analysis_delete_delta)
            if project_id == "": continue
            if project_id not in good_records_dict:
                good_records_dict[project_id] = {"delivery_start_time":delivery_start_time, "ID":[]}
            good_records_dict[project_id]["ID"].append(id)
            delivery_start_time_old = good_records_dict[project_id]["delivery_start_time"]
            time_delta = self.time_handle.delta_days(delivery_start_time-delivery_start_time_old)
            #如果该条记录时间比上条记录时间更靠近现在，则将delivery_start_time时间更新为该条记录的启动交付时间，否则保留原启动交付时间。
            #即同一个项目符合条件的记录中，保留最新的启动交付时间
            if time_delta > 0 :
                good_records_dict[project_id]["delivery_start_time"] = delivery_start_time
            else:
                continue

         #如果项目记录有多个，以最近一次的为准
        return good_records_dict
    def get_all_records(self):
        '''
        功能：获取数据库中 交付是否完成=是,是否删除=否 的项目列表
        return:所有项目列表，list类型[(),(),()]
        '''
        select_conditions = [ ("delivery_bool",1), ("local_delete_bool",0) ]
        tb_name = self.config_dict["table"]
        all_records = self.local_sql.select( table_name=tb_name, select_conditions=select_conditions )
        if len(all_records ) == 0 :
            my_log.info( "没有需要删除的项目记录" )
            return []
        else:
            return all_records
    
    def get_good_record( self, tmp ):
        '''
        功能：根据提供的纪录，判断是否满足删除要求
        过程：
            如果删除类型是超期删除，则统一按照过滤项目的删除时间判断（因为如果是分析项目，超过21天的时候提醒插入到删除库后，按照30天判断，则会在51天的时候删除，讨论后认为不合理）
            如果删除类型是交付删除：
                如果是分析项目，则按照分析项目的删除时间判断，满足条件，则返回项目编号，启动交付时间和待删除数据库中的id
                如果是非分析项目，则按照非分析项目的删除时间判断，满足条件，则返回项目编号，启动交付时间和待删除数据库中的id
        返回值:
            project_id, delivery_start_time, id (分别为子项目编号，启动交付时间和待删除数据库中的id)
        '''
        id = tmp[0]
        self.project_id = tmp[1]
        self.delivery_start_time = tmp[2] #启动交付时间
        self.analysis_bool = tmp[4] #是否为分析项目 0-否，1-是
        self.delete_type = tmp[5]
        self.delivery_end_time = tmp[9] 
        self.local_delete_bool = tmp[10] #本地是否删除
        time_handle = TimeFormat(self.delivery_end_time, now_time )
        now_time = time_handle.get_now_time(self.date_format)
        self.delta_days = time_handle.get_delta_days()
        #如果是超期删除，则将该项目的analysis_bool是否为分析项目的值改成0，即非分析项目，后面按照非分析项目的删除规则执行判断
        if self.delete_type == "超期删除":
            self.analysis_bool = 0 
        if self.analysis_bool == 1 :
            if self.delta_days > self.filter_delete_delta :
                return self.project_id, self.delivery_start_time, id
            else:
                return "","",""
        else :
            if self.delta_days > self.analysis_delete_delta :
                return self.project_id, self.delivery_start_time, id
            else:
                return "","",""
                      
class DeleteFQ():
    def __init__( self, project_dict, config_dict, lims_db, local_sql, test_mode ):
        self.project_dict = project_dict
        self.config_dict = config_dict
        self.lims_db = lims_db
        self.local_sql = local_sql
        self.test_mode = test_mode
        self.time_handle = TimeFormat()
        self.retry_times = 3
        self.fail_list = []
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
        for project_id in self.project_dict :
            delivery_start_time = self.project_dict[project_id]['delivery_start_time']
            id_list = self.project_dict[project_id]['ID']
            project_filter_path = self.get_filter_dir( project_id ) #整个目录，即Analysis_1/../
            filter_path_glob_list = glob.glob( project_filter_path+'/*' ) #获取项目过滤路径下的所有路径和文件
            filter_path_list = self.filter_glob(filter_path_glob_list) #将文件和Filter_Result过滤掉,只保留Anlysis*开头的目录
            if len(filter_path_list) ==0 : continue
            #一个项目一个删除列表
            delete_list = []
            for filter_analysis_dir in filter_path_list:
                self.get_delete_list( filter_analysis_dir, project_id, delivery_start_time, delete_list )
            if delete_list :
                self.run_delete( delete_list, project_id, delivery_start_time )
            else:
                #为什么有待删除记录，但是没有要删除的文件呢？ 需要查看
                my_log.error( "项目 {0} 没有需要删除的文件".format(project_id ))
        return self.fail_list
    def get_delete_list(self, filter_dir, project_id, delivery_start_time, delete_list ):
        '''
        功能：获取删除列表，并将其删除
        过程：
            第一步：获取过滤路径创建时间是否在启动交付时间之前，如果是，则继续判断
            第二步：通过get_qc_dir和get_chip_dir分别获取该Analysis路径对应芯片的下机数据路径和Filter_Result目录下的CHIP/CHIP_FC
            第三步：如果以上两个路径不为空，则添加到待删除列表中
        参数：
            filter_dir:过滤Analysis_1目录
            project_id:子项目编号
        返回值：无
        '''
        retry_times = 3
        dir_create_time = self.time_handle.dit_create_time(filter_dir) #获取目录的创建时间
            
        time_delta = self.time_handle.time_delta(dir_create_time, delivery_start_time) #获取启动交付时间-目录创建时间
        #如果time_delta>=0 说明目录创建时间不迟于数据启动交付时间,需要删除,如果小于0，说明启动交付时间早于目录创建时间，不需要删除
        if time_delta > 0 :
            delete_list.append( filter_dir )
            qc_dir = self.get_qc_dir( project_id, filter_dir)
            chip_dir = self.get_chip_dir( filter_dir )
            if qc_dir != '' : 
                delete_list.append( qc_dir )
            if chip_dir != '' :
                delete_list.append( chip_dir )
        return delete_list
    def run_delete( self, delete_list, id_list ):
        '''
        功能：
            对删除列表执行删除操作。
            如果删除成功，则将删除数据库的记录更新[("local_delete_bool",1),('delete_time',now_time)]，如果更新失败，输出到.e文件中
            如果删除失败，再尝试3次，如果仍然失败，则输出到.e文件中
        注意事项：如果是测试模式（self.test_mode=True），则不执行删除操作且不更新数据库，只输出删除命令
        参数：
            delete_list:删除列表，list类型，格式为[,,]
            id_list：待删除项目在删除数据库中的id信息，用于更新删除数据库中的记录
            self.retry_times:删除命令失败的重试次数
            self.test_mode：是否为测试模式
        '''
        delete_status = []
        for delete_dir in delete_list:
            retry_times = self.retry_times
            delete_cmd = 'rm -rf {0}'.format(delete_dir)
            run_status = cmd_run(delete_cmd, self.test_mode)
            #如果删除命令失败，重试3次
            while( run_status==False & retry_times > 0 ) :
                run_status = cmd_run(delete_cmd, self.test_mode )
                retry_times -= 1
            if run_status == False:
                delete_status.append( False )
                my_log.error("{0} Failed".format( delete_cmd))
            else:
                delete_status.append( True )
        if set( delete_status ) :
                #如果是测试模式不更新删除数据库状态
            if not self.test_mode : 
                now_time = self.time_handle.now_time()
                update_info = [("local_delete_bool",1),('delete_time',now_time)]
                for each_id in id_list:
                    update_condition = [("ID", each_id)]
                    try:
                        self.local_sql.update(self.config_dict['table'],value_list=update_info, conditions=update_condition )
                    except:
                        my_log.error("update {0} in {1} Failed".format( each_id, self.config_dict['table'] ))
        else:
            self.fail_list.apeend(','.join(id_list))
    def get_filter_dir(self, project_id):
        '''
        功能：根据项目编号从tb_filter_task表中获取项目过滤路径
        注意事项：默认tb_fitler_task中的过滤路径都是一样的，所以只获取tb_filter_task表中的第一条记录
        参数：
            project_id:项目编号
        返回值：过滤路径（为tb_filter_task中的 analysis_path/../路径 
        '''
        select_conditions = [('merge_project_code',project_id)]
        filter_records = self.lims_db.select( table_name='tb_filter_task', col_list=["analysis_path"], conditions=select_conditions) 
        if len(filter_records) == 0:
            my_log.error("{0} not in tb_filter_task")
        else:
            analysis_path = filter_records[0][0] #获取第一个记录的分析路径
            filter_path = os.path.abspath("{0}/../".format(analysis_path))
            return filter_path

    def filter_glob(self, filter_path_glob_list ):
        '''
        功能：将过滤路径下的文件筛选掉，只保留目录，文件可能是临时文件或者是一些记录文件
        参数：
            filter_path_glob_list：过滤路径下的文件/目录列表以及Filter_Result目录
        return:path_list(只有目录)
        '''  
        path_list = []
        for path in filter_path_glob_list:
            if os.path.isdir(path):
                if os.path.basename(path).startswith('Analysis'):
                    path_list.append( path )
        return path_list
    
    def get_qc_dir(self, project_id, filter_path ):
        '''
        功能：根据过滤的Analysis_1_FC路径获取对应的芯片号，并获取该批次的拆分目录。如果目录不是以Analysis开头,则返回空
        参数：
            filter_path:过滤的Analysis_1_FC路径
        return:
            qc_dir:拆分目录
        '''
        fc_no = self.get_fc( filter_path )
        select_conditions = [('merge_project_code',project_id), ('fc_no',fc_no),("deleted","NO"), ("filter_status","FILTER_FINISHED") ]
        filter_qc_dir_records = self.lims_db.select( table_name='tb_filter_task', col_list=["raw_dir"], conditions=select_conditions) 
        if len(filter_qc_dir_records) == 0:
            my_log.error("{0} {1}not in tb_filter_task".format( project_id, fc_no))
            qc_dir = ""
        else:
            qc_dir = filter_qc_dir_records[0][0] #获取第一个记录的分析路径
        return qc_dir 
    def get_chip_dir( self, filter_path ):
        '''
        功能：获取过滤Filter_Result目录下的CHIP目录下的芯片数据，避免Analysis_1_FC删除后，该目录数据仍然存在，数据过滤时会报错
        参数：
            filter_path:过滤的Analysis_1_FC路径
        return:
            fc_chip_dir:CHIP目录
        '''
        fc_no = self.get_fc( filter_path )
        fc_chip_dir = '{0}/../Filter_Result/CHIP/CHIP_{1}/'.format( filter_path, fc_no )
        if os.path.abspath( fc_chip_dir ): 
            return fc_chip_dir
        else:
            return ""
    def get_fc( self, filter_path ):
        filter_analysis_path = os.path.abspath(filter_path)
        filter_analysis_path_end = filter_analysis_path.split('/')[-1]
        if filter_analysis_path_end.startswith('Analysis'):
            fc_no = filter_analysis_path_end.split('_')[-1]
            return fc_no
        return ""


def cmd_run( cmd, test_mode ):
    '''
    功能：执行os.system命令的操作，如果是test_mode模式，则只打印命令，返回True
    参数：
        cmd:执行命令
        test_mode:是否为test_mode模式，True为test_mode模式，False为正式模式;如果为True,不执行命令，返回True;如果为False，执行命令
    返回值：bool,True or False
    '''
    if test_mode:
        my_log.info("Test:{0}".format(cmd))
        return True
    else :
        if os.system(cmd) == 0:
            my_log.info("{0} Success".format( cmd))
            return True
        else:
            return False
def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument('-c','--pipe_config',help='config_file',dest='pipe_config',default = '{0}/../config/config.txt'.format(bindir))
    args=parser.parse_args()

    print(args.pipe_config)
    date_format= '%Y-%m-%d %H:%M:%S'
    time_handle = TimeFormat( date_format )
    now_time = time_handle.now_time()
    #读取config文件以及数据库连接
    config_dict = DOMysql.read_config( args.pipe_config )
    lims_db = DOMysql.SQL( args.pipe_config, 'lims' )
    local_sql = DOMysql.SQL( args.pipe_config )
    
    select_record = SelectItem(local_sql, config_dict, date_format)
    project_record = select_record.main() #{project_id:{"delivery_start_time":"2024-05-17 09:43:56", "id":[1,2,3]}}
    delete_handle = DeleteFQ( project_record, config_dict, lims_db, local_sql, args.test )
    fail_list = delete_handle.main()
    #删除失败的机器人发消息到群里
    #当前的机器人是 数据删除讨论群 里的机器人
    robot_url = "https://qyapi.weixin.qq.com/cgi-bin/webhook/send?key=ffa61817-1dbe-4d94-981b-b3572abfe470"
    title = "删除目录失败的id"
    label = "id"
    content = [[i] for i in fail_list]
    robot.Robots(title=title,label_list=label,content_list=content, url=robot_url)


if __name__ == '__main__':
    main()


