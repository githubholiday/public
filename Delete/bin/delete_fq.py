#-*- encoding:utf-8 -*-
'''
【功能】
根据删除数据库中的信息获取符合条件的项目数据（过滤路径（Analysis和CHIP）和拆分路径），并输出到指定文件中
获取配置文件中的outdir路径，并在路径下输出：
1）日期_filter_rm.sh：过滤路径下符合删除条件的目录命令
2）日期_sci_qc_rm.sh：过滤路径下符合删除条件的目录对应芯片的拆分路径的删除命令
3) 日期_notice_info.txt：需要引起重视，去查看的文件。如项目的完成交付时间满足7天，但是过滤路径下没有符合创建时间早于启动交付时间的目录，需要排查原因，是否要升级脚本
4) 日期_record.txt: 本地删除的项目、过滤路径，启动交付时间，目录创建时间的记录信息

【规则】
1. 筛选交付完成&未删除的项目条目
2. 如果交付完成时间超过指定时间（过滤7天，非分析30天），则认为可以删除
3. 如果项目过滤路径下的目录创建时间早于启动交付时间，则将该项目的该过滤目录、Filter_Result/CHIP、拆分目录均删除

【操作步骤】
1. 筛选符合条件的项目列表
交付是否完成=是
是否删除=否

2. 计算1中所有项目列表的 交付完成时间 与 当前时间 的差值天数 delta_days
如果是删除类型=超期删除的项目，则按照非分析项目的判断条件执行
1)如果是分析项目，delta_days >30天(analysis_delete_delta),存入待删除字典
2)如果是非分析项目,delta_days>7天(配置文件中的filter_delete_delta),存入待删除字典
=>获取待删除字典 good_records_dict[project_id] = {"delivery_start_time":delivery_start_time, "ID":[1,2,3,4]} 
=>如果同一个项目有多条记录，则启动交付时间获取的最近的一次，即启动交付时间为最新一次的未删除记录的启动交付时间
=>如果待删除列表为空，则屏幕输出无记录

3. 根据2中的待删除列表，获取项目的过滤路径，判断过滤路径下的目录是否有在启动交付时间前创建的，如果满足条件，则为可删除
1）通过2中的待删除列表的项目编号，从tb_filter_task列表中，获取到该项目的过滤路径 （tb_filter_task.filter_path/../）
2）glob过滤路径下的所有目录，包括Analysis_1和Filter_Result，只筛选Analysis*开头的目录，用于下一步分析
3）获取每个目录的创建时间
如果 目录创建时间早于启动交付时间，(函数：get_delete_list -启动交付时间-创建时间>0）
3.1）将过滤路径删除命令存入到日期_filter_rm.sh文件中
3.2）获取其对应的芯片号，并从tb_filter_task中获取对应的拆分数据路径和Filter_Result/CHIP/CHIP_FC路径，分别存入到 日期_sci_qc_rm.sh 和 日期_filter_rm.sh文件中
如果该项目的所有过滤目录时间都晚于启动交付时间，则将项目信息输出到 日期_notice_info.txt 中，需要查看原因

4）根据good_records_dict中的id_list逐条更新删除数据库中该项目的删除状态local_delete_bool=1,delete_time=当前时间(年-月-日) 


【参数说明】
-c,--pipe_config:配置文件，使用到配置文件中的信息有 数据库信息,删除时间中的filter_delete_time
如果是测试版本，则使用config/config_test.txt文件

【注意事项】
1. 删除项目时，如果同一个项目有多条未删除记录，则以最新的启动交付时间为准
2. 输出文件中的过滤和拆分目录需要分别使用对应的账号进行删除，并check是否有删除失败的，如果有删除失败的，需要将其重新删除，直到成功为止（因为删除数据库已经更新为删除成功了）
3. 过滤项目重过滤的时候，将Analysis目录删除重新创建，使得目录的创建时间有更新
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
		#create_stamptime = os.path.getmtime(dir_path) #时间戳 163467344
		#create_localtime= time.localtime(create_stamptime)
		#create_formatime = time.strftime(self.date_format ,create_localtime) #格式化为时间 2024-04-14
		#2024-7-10修改，上面获取的目录创建时间是错误的
		#create_formatime = datetime.fromtimestamp(os.stat(dir_path).st_mtime).strftime(self.date_format)
		create_formatime = datetime.fromtimestamp(os.stat(dir_path).st_atime).strftime(self.date_format)

		return create_formatime
	def time_stamp_transfer(self, time_stamp, date_format ='%Y-%m-%d %H:%M:%S'):
		'''
		功能：将时间戳转化为时间格式
		参数:
			time_stamp :时间戳，如1714685197412
			date_format:时间格式，默认为 '%Y-%m-%d %H:%M:%S'
		'''
		timestamp = time.localtime(time_stamp/1000) #此处是否除以1000，要看你的时间戳单位，此处保持跟上方一致
		time_format = time.strftime(date_format,timestamp)
		return time_format

class SelectItem():
	def __init__(self, local_sql, config_dict, date_format = '%Y-%m-%d %Y-%m-%d'):
		self.local_sql = local_sql
		self.date_format = date_format
		self.config_dict = config_dict
		self.filter_delete_delta = int(config_dict['filter_delete_delta'])
		self.analysis_delete_delta = int(config_dict['analysis_delete_delta'])
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
			project_id, delivery_start_time, id  = self.get_good_record( tmp )
			if project_id == "": continue
			if project_id not in good_records_dict:
				good_records_dict[project_id] = {"delivery_start_time":delivery_start_time, "ID":[]}
			good_records_dict[project_id]["ID"].append(id)
			delivery_start_time_old = good_records_dict[project_id]["delivery_start_time"]
			time_delta_seconds = self.time_handle.get_delta_seconds(delivery_start_time_old,delivery_start_time)
			#如果该条记录时间比上条记录时间更靠近现在，则将delivery_start_time时间更新为该条记录的启动交付时间，否则保留原启动交付时间。
			#即同一个项目符合条件的记录中，保留最新的启动交付时间
			if time_delta_seconds > 0 :
				good_records_dict[project_id]["delivery_start_time"] = delivery_start_time
			else:
				continue
		if len(good_records_dict) == 0 :
			my_log.info("无记录-数据删除数据库中没有交付完成时间的超过指定时间的项目记录,过滤项目-{self.filter_delete_delta},分析项目-{self.analysis_delete_delta}".format( self=self))
		 #如果项目记录有多个，以最近一次的为准
		return good_records_dict
	def get_all_records(self):
		'''
		功能：获取数据库中 交付是否完成=是,是否删除=否 的项目列表
		return:所有项目列表，list类型[(),(),()]
		'''
		select_conditions = [ ("delivery_bool",1), ("local_delete_bool",0) ]
		tb_name = self.config_dict["table"]
		all_records = self.local_sql.select( table_name=tb_name, conditions=select_conditions )
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
				如果是分析项目，则按照分析项目的删除时间判断，现在-交付完成时间满足条件，则返回项目编号，启动交付时间和待删除数据库中的id
				如果是非分析项目，则按照非分析项目的删除时间判断，现在-交付完成时间满足条件，则返回项目编号，启动交付时间和待删除数据库中的id
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
		now_time = self.time_handle.now_time()
		#获取当前时间-交付结束时间的天数
		self.delta_days = self.time_handle.get_delta_days(self.delivery_end_time,now_time)# now_time)
		#如果是超期删除，则将该项目的analysis_bool是否为分析项目的值改成0，即非分析项目，后面按照非分析项目的删除规则执行判断
		if self.delete_type == "超期删除":
			self.analysis_bool = 0 
		#如果是分析项目
		if self.analysis_bool == 1 :
			if self.delta_days > self.analysis_delete_delta :
				return self.project_id, self.delivery_start_time, id
			else:
				return "","",""
		#如果是过滤项目
		else :
			if self.delta_days > self.filter_delete_delta :
				return self.project_id, self.delivery_start_time, id
			else:
				return "","",""
					  
class DeleteFQ():
	def __init__( self, project_dict, config_dict, lims_db, local_sql, filter_sh, sci_qc_sh, notice_file, record_file ):
		self.project_dict = project_dict
		self.config_dict = config_dict
		self.lims_db = lims_db
		self.local_sql = local_sql
		self.filter_rm_sh = filter_sh
		self.sci_qc_rm_sh = sci_qc_sh
		self.notice_file = notice_file
		self.record_file = record_file
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
			#获取tb_filter_task中在本地仍然存在的过滤路径以及过滤开始时间
			project_filter_path_dict = self.get_filter_dir( project_id )
			if not project_filter_path_dict :
				self.notice_file.write("*error4*{0} 信息在数据删除数据库中并符合删除逻辑,但是不在tb_filter_task数据库中,.可能是本地路径已经不存在了,需要排查\n".format(project_id))
				continue
			#一个项目一个删除列表,如果有一个目录符合条件能删除，就可以更新删除数据库的删除状态；如果没有目录符合条件，则不更新删除数据库的删除状态，并输出到.e文件中
			delete_bool_list = []
			for filter_analysis_dir in project_filter_path_dict:
				filter_path_create_time = project_filter_path_dict[filter_analysis_dir]
				delete_bool = self.get_delete_list( filter_analysis_dir, project_id, delivery_start_time, filter_path_create_time )
				delete_bool_list.append(delete_bool)
			if True in set(delete_bool_list) :
				self.update_db( id_list  )
			else:
				self.notice_file.write("*error5* {0} 过滤目录下没有符合条件的删除记录:均不符合目录的开始过滤时间早于启动交付时间\n".format(project_id))
				for filter_analysis_dir in project_filter_path_dict:
					filter_path_create_time = project_filter_path_dict[filter_analysis_dir]
					self.notice_file.write("*error5* {0} {1} 过滤开始时间{2} 晚于启动交付时间 {3}\n".format(project_id,filter_analysis_dir,  filter_path_create_time, delivery_start_time))


	def get_delete_list(self, filter_dir, project_id, delivery_start_time, dir_create_time ):
		'''
		功能：获取删除列表，并将其删除
		过程：
			第一步：获取过滤路径创建时间是否在启动交付时间之前，如果是，则继续判断
			第二步：通过get_qc_dir和get_chip_dir分别获取该Analysis路径对应芯片的下机数据路径和Filter_Result目录下的CHIP/CHIP_FC
			第三步：将过滤路径、CHIP路径和芯片路径写入到filter_rm_sh文件和sci_qc_rm_sh文件
		参数：
			filter_dir:过滤Analysis_1目录
			project_id:子项目编号
			delivery_start_time:启动交付时间
			dir_create_time: 过滤路径(Analysis_1_)创建时间
		返回值：无
		'''
		retry_times = 3
		#dir_create_time = self.time_handle.dir_create_time(filter_dir) #获取目录的创建时间
		time_delta_seconds = self.time_handle.get_delta_seconds(dir_create_time,delivery_start_time )
		#如果time_delta>=0 说明目录创建时间不迟于数据启动交付时间,需要删除,如果小于0，说明启动交付时间早于目录创建时间，不需要删除
		if time_delta_seconds > 0 :
			#my_log.info("{0} 满足删除条件，创建时间早于启动交付时间".format( filter_dir ))
			self.filter_rm_sh.write( 'rm -rf ' + filter_dir+'\n' )
			qc_dir = self.get_qc_dir( project_id, filter_dir)
			chip_dir = self.get_chip_dir( filter_dir )
			self.record_file.write("{project_id}\t{filter_dir}\t{delivery_start_time}\t{dir_create_time}\n".format(project_id=project_id, filter_dir=filter_dir, delivery_start_time=delivery_start_time, dir_create_time=dir_create_time ))
			if chip_dir :
				self.filter_rm_sh.write( 'rm -rf ' + chip_dir+'\n' )
			else:
				self.notice_file.write("*error1* {0} 过滤目录下没有对应的CHIP目录:{1}\n".format(project_id,filter_dir ))
			if qc_dir :
				self.sci_qc_rm_sh.write( 'rm -rf ' +qc_dir + '\n' )
			else:
				self.notice_file.write("*error2* {0} 没有对应的拆分数据路径:{1}\n".format(project_id,filter_dir ))
			return True
		else:
			return False

	def update_db( self, id_list ):
		'''
		功能：
			根据id将删除数据库中的条目状态更改
		参数：
			id_list：删除数据库中的id列表
		'''
		for each_id in id_list:
			update_condition = [("ID", each_id)]
			now_time = self.time_handle.now_time()
			update_info = [("local_delete_bool",1),('delete_time',now_time)]

			try:
				self.local_sql.update(self.config_dict['table'],value_list=update_info, conditions=update_condition )
				#my_log.info("Update {0} in {1} Successfully".format( each_id, self.config_dict['table'] ))
			except:
				my_log.error("update {0} in {1} Failed".format( each_id, self.config_dict['table'] ))

	def get_filter_dir(self, project_id):
		'''
		功能：根据项目编号从tb_filter_task表中获取项目过滤路径
		注意事项：默认tb_fitler_task中的过滤路径都是一样的，所以只获取tb_filter_task表中的第一条记录
		参数：
			project_id:项目编号
		返回值：过滤路径（为tb_filter_task中的 analysis_path/../路径 
		'''
		select_conditions = [('merge_project_code',project_id),("deleted","NO")]
		filter_records = self.lims_db.select( table_name='tb_filter_task', col_list=["analysis_path","filter_start_time"], conditions=select_conditions) 
		filter_path_dict = {}
		if len(filter_records) == 0:
			my_log.error("{0} not in tb_filter_task".format(project_id))
			return ""
		for filter_info in filter_records:
			if filter_info[0] == "":continue
			filter_path = os.path.abspath(filter_info[0])
			##2024-8-14 减少错误
			if len(filter_path.split('/')) <= 10: 
				self.notice_file.write("{0} 过滤路径格式不满足要求".format(filter_path))
				continue
			filter_start_time = int(filter_info[1])
			if os.path.exists( filter_path ) :
				if filter_path not in filter_path_dict:
					filter_path_dict[filter_path] = self.time_handle.time_stamp_transfer( filter_start_time )
				else:
					self.notice_file.write("*error3* {0} {1} 在tb_filter_task中有多个符合条件的记录，请核实\n".format( project_id, filter_path ))
		return filter_path_dict
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
			my_log.error("{0} {1} not in tb_filter_task".format( project_id, fc_no))
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
		fc_chip_dir = os.path.abspath('{0}/../Filter_Result/CHIP/CHIP_{1}/'.format( filter_path, fc_no ))
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

	my_log.info("配置文件:"+args.pipe_config)
	#读取config文件以及数据库连接
	my_log.info("连接lims数据库")
	lims_db = DOMysql.SQL( args.pipe_config, db_type='lims' )
	my_log.info("连接数据删除数据库")
	local_sql = DOMysql.SQL( args.pipe_config )
	date_format = local_sql.config_dic['TimeFormat']
	#
	time_handle = TimeFormat( "%Y-%m-%d-%H-%M-%S" )
	now_time = time_handle.now_time()
	now_date = TimeFormat( "%Y-%m-%d" ).now_time()

	# 
	pre_outdir = local_sql.config_dic['outdir']
	outdir = '{0}/{1}'.format( pre_outdir, now_date )
	if not os.path.exists( outdir ):
		os.makedirs( outdir )
	filter_rm_sh = '{0}/{1}_filter_rm.sh'.format( outdir, now_time )
	sci_qc_rm_sh = '{0}/{1}_sci_qc_rm.sh'.format( outdir, now_time )
	notice_file = '{0}/{1}_notice_info.txt'.format( outdir, now_time )
	record_file = '{0}/{1}_record.txt'.format( outdir, now_time )
	filter_rm_sh_handle = open( filter_rm_sh , 'w')
	sci_qc_rm_sh_handle = open( sci_qc_rm_sh, 'w')
	notice_file_handle = open( notice_file, 'w')
	record_file_handle = open( record_file, 'w')
	record_file_handle.write('\t'.join(['项目编号','过滤路径','启动交付时间','目录创建时间']) + '\n')

	#筛选交付完成时间符合条件的（即交付完成距离当前时间超过 过滤和分析项目删除周期，一般为7和30天)
	select_record = SelectItem(local_sql, local_sql.config_dic, date_format)
	project_record = select_record.main() #{project_id:{"delivery_start_time":"2024-05-17 09:43:56", "id":[1,2,3]}}
	my_log.info("获取交付完成时间满足删除周期的项目列表")
	#对上面筛选的项目列表 判断过滤路径是否需要删除（过滤路径创建时间是否早于启动交付时间），如果需要删除，则分别输出到 filter_rm_sh_handle, sci_qc_rm_sh_handle
	#notice_file_handle:是一些需要关注的信息，例如项目路径下没有符合要求的项目，需要查找原因，其余内容均屏幕输出了
	delete_handle = DeleteFQ( project_record, local_sql.config_dic, lims_db, local_sql, filter_rm_sh_handle, sci_qc_rm_sh_handle, notice_file_handle, record_file_handle )
	fail_list = delete_handle.main()

	filter_rm_sh_handle.close()
	sci_qc_rm_sh_handle.close()
	notice_file_handle.close()
	record_file_handle.close()
	my_log.info("过滤路径删除命令文件: {0}".format(filter_rm_sh))
	my_log.info("拆分路径删除命令文件: {0}".format( sci_qc_rm_sh))
	my_log.info("***需要查看是否异常的文件: {0}".format( notice_file))
	my_log.info("删除记录表:{0}".format( record_file ))
	#删除失败的机器人发消息到群里
	#当前的机器人是 数据删除讨论群 里的机器人
	
	robot_url = pre_outdir = local_sql.config_dic['robot_url'] 
	title = "{0} 删除过滤目录脚本已经执行完毕,请删除拆分和过滤数据,并查看异常列表".format( now_time )
	label = ['日期','输出路径']
	content = [now_time, outdir ]
	robot.Robots(title=title,label_list=label,content_list=content, url=robot_url)
	filter_rm_sh_handle.close()
	sci_qc_rm_sh_handle.close()


if __name__ == '__main__':
	main()



