#! /usr/bin/env python3
# -*- coding: utf-8 -*-  
'''
说明：
	任务状态有八种：
	* waiting : 已生成文件，但没有投递
	* qsub：已投递，尚未运行
	* running： 运行中
	* hold： hold住了
	* done： 在之前就跑完了，这次不用跑
	* finished： 这次跑完了
	* break： 断掉了
	* fail： 怎么投都错了，所以放弃了
程序自动退出的情况：
	* 某个任务反复投递后，都没有运行成功，报错为 ：Exit : {0} is not finish
	* 某个命令(不包括du和qstat -j)反复运行出错，导致程序退出： 报错为：Exit : Command {0} is error
	* qsub反复投递不上，报错为：Exit :  qsub failed {0}
	* 无法创建目录，报错为: Exit : cannot mkdir -p {0}
	* 数据库查询失败，报错为：Exit ： cannot access to db
程序发现异常情况后，如何排查：
1. 检查所有的.e和.o文件是否有报错出现
	1.1 如果.e中有报错，则对应的去查错
	1.2 如果.e中无报错，.o中也无 This-Work-is-Completed! ，则如果多次投递都是异常退出，则可能是 节点中断造成的,可以去work.sh.log中查看任务运行的节点，是否有一定规定。
2. 检查 work.sh.jobid.log文件中是否有反复无法投递的任务
3. 检查work.sh.log文件中是否有退出的信息,该文件中信息很多，直接搜索关键词： WARNING ，Error, Exit等，INFO类基本不用关注。
4. 检查最顶层的nohup.out 是否有python退出的报错信息。


更新说明：
2015年12月22日：
1. 使用subprocess中的wait来阻塞进程，防止qstat无法响应的状态，并使用poll检查返回值
2. 修复了failed状态又重新被修改成break的状态。
3. 增加了状态返回值，如果全部运行完，则为exit(0),否则为exit(1)
2016年3月30日
1. 修改了du 无响应的返回值，返回为[]， 不退出程序
2. 降低了qhost的频率
3. 降低了qstat -j 的频率
2016年7月14日：
1. 修复了任务数太多，无法投递上去的bug。如果太多，会等待可以投递的时候再投递。
2016年7月29日：
1. 修复了qstat过多，对系统造成的负载。增加了qstat无响应，等待5min，重复10次；并且每次qstat后，sleep 5s，减少负载。
2016年8月1日
1. 修复了任务运行过快，导致了qstat时存在，但qstat -j却不存在的bug
2016年8月18日：
1. 增加了popen的timeout参数，以对抗相应超时的时候，可以自动结束进程。du的相应时间设为300s；其他的是60s；
2016年11月29日：
1. 修复了.o文件可能是乱码，导致检查出错的情况。
2. 优化了debug_log日志，可以更方便的查错
2016年12月26日:
1. 增加了-js/--jobidstart参数，来指定job id count的起始数，默认是0,以适配流程中的-a参数
'''
import argparse
import sys
import os
import re
import time
import glob
import logging
import subprocess
import random
import threading
import sqlite3
from datetime import datetime
import drmaa
bindir = os.path.abspath(os.path.dirname(__file__))

MEMORY_INCREASE_RATIO = 1.5



#qhost_counter = 0 
__author__='Liu Tao'
__mail__= 'taoliu@annoroad.com'

pat1=re.compile('^\s+$')
def mylogger(script):
	logger = logging.getLogger('mylogger')
	logger.setLevel(logging.DEBUG)

	fh = logging.FileHandler('{0}.log'.format(script))
	fh.setLevel(logging.DEBUG)

	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	return logger

class Command():
	def __init__(self, cmd):
		self.cmd = cmd
		self.process = None
		self.content = '' 
	def run(self, timeout): ## https://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
		def target():
			#print ('Thread started')
			self.process = subprocess.Popen(self.cmd, shell=True , stdout = subprocess.PIPE)
			#self.process.wait()
			self.content = self.process.communicate()[0]
			#self.content = [ i for i in self.process.stdout]
			#print ('Thread finished')

		thread = threading.Thread(target=target)
		thread.start()

		thread.join(timeout)
		if thread.is_alive():
			print ('Terminating process')
			self.process.terminate()
			thread.join()
		if self.process == None :
			# in case for BlockingIOError: [Errno 11] Resource temporarily unavailable in subprocess.Popen 
			# make self.process is None also
			return (1)
		else :
			return ( self.process.returncode)

class Database():
	def __init__(self, path ):
		self.database = path
		self.conn = sqlite3.connect( self.database )
	def check_table_exists(self):
		cmd = 'select name from sqlite_master where type = "table";'
		cursor = self.execute(cmd)
		tables = [ j for i in cursor.fetchall() for j in i]
		r_value = True
		debug_log.info('table: {0}'.format(tables))
		for i in [ 'PARAMETER' , 'SCRIPT' , 'JOB' ]:
			if len(tables) ==  0 or  not i in tables:
				self.create_database( i )
				r_value = False
		if len(self.query_all('SCRIPT')) == 0 : 
			r_value = False 
		return r_value
	def create_database(self , type ):
		cmd = ''
		if type == 'PARAMETER':
			cmd = '''CREATE TABLE PARAMETER
				(ID        TEXT    PRIMARY KEY,
				MAXJOBS    TEXT    ,
				INPUT   TEXT ,
				QUOTA     TEXT      ) ; '''
		elif type == 'SCRIPT':
			cmd = '''CREATE TABLE SCRIPT 
				(ID       TEXT      PRIMARY KEY  ,
				 NAME     TEXT     NOT NULL ,
				 SCRIPT   TEXT ,
				 QUEUE    TEXT , 
				 P        TEXT , 
				 VF       TEXT ,
				 HVMEM    TEXT , 
				 OTHER    TEXT , 
				 OLDP        TEXT , 
				 OLDVF       TEXT ,
				 OLDHVMEM    TEXT , 
				 MAXCYCLE TEXT , 
				 JOBID    INT ,
				 QSUBTIME INT , 
				 FINISH   INT ,
				 STATE    TEXT) ; '''
		elif type == 'JOB':
			cmd = '''CREATE TABLE JOB 
				(JOBID     INT       PRIMARY KEY  , 
				 NAME      TEXT      NOT NULL , 
				 SCRIPT    TEXT      NOT NULL , 
				 STATE     TEXT      NOT NULL , 
				 P         TEXT      NOT NULL , 
				 VF        TEXT      NOT NULL , 
				 HVMEM     TEXT      NOT NULL , 
				 OTHER     TEXT      NOT NULL , 
				 NODE      TEXT      NOT NULL , 
				 UCPU      INT       NOT NULL ,
				 UMEM      TEXT      NOT NULL , 
				 UIO       TEXT      NOT NULL , 
				 UMAXVMEM   TEXT   ,
				 UTHREAD   TEXT   ,
				 START     TEXT   ,
				 END       TEXT   ,
				 FINISH    TEXT   ,
				 KILL      TEXT   ,
				 DIED      TEXT ); '''
		self.execute(cmd)
		debug_log.info(cmd)
		self.commit()
	def query(self , table_name , key , a_value , cols):
		cmd = 'SELECT {3} FROM {0} WHERE {1} = "{2}"; '.format(table_name , key , a_value , ",".join(cols))
		debug_log.info(cmd)
		cursor = self.execute(cmd)
		tt = cursor.fetchall()
		return tt
	def query_all(self , table_name):
		cmd = 'SELECT * FROM {0};'.format(table_name)
		debug_log.info(cmd)
		cursor = self.execute(cmd)
		tt = cursor.fetchall()
		return tt
	def insert(self , table_name , name_list , value_list):
		names = ",".join(name_list) 
		values = ",".join(['"' + str(i) + '"' for i in  value_list])
		cmd = 'INSERT INTO {0} ({1}) VALUES ({2});'.format(table_name , names , values)
		debug_log.info(cmd)
		self.execute(cmd)
		self.commit()
	def delete(self , table_name , key , a_value):
		cmd = 'DELETE FROM {0} WHERE {1} = {2} ; '.format(table_name , key , a_value)
		debug_log.info(cmd)
		self.execute(cmd)
		self.commit()
	def delete_all_content(self , table_name ):
		cmd = 'DELETE  FROM {0} ;  '.format(table_name )
		debug_log.info(cmd)
		self.execute(cmd)
		self.commit()
	def update(self , table_name , record , key , a_value):
		cmd = ''
		for i , j in record.items():
			cmd = 'UPDATE {0} set {1} = "{2}" WHERE {3} = "{4}";'.format(table_name , i, j , key , a_value )
			debug_log.info(cmd)
			self.execute(cmd)
			self.commit()
	def commit(self):
		self.conn.commit()
	def execute(self , cmd):
		count = 10
		while count :
			try :
				return self.conn.execute(cmd)
			except sqlite3.OperationalError : 
				count -= 1
				time.sleep(10)
		else:
			sys.exit('sqlite3.OperationalError:  unable to open database file')
	def close(self):
		self.conn.close()

'''
To solve one problem in commander:
	error: commlib error: got select error (Connection refused)
	Unable to run job: unable to send message to qmaster using port 6444 on host "bjannoroad.loc
	al": got send error.
	Exiting.

20180209
'''
def Submit(FileName, para):
	try:
		debug_log.info('qsub prepared ' )
		with drmaa.Session() as s : 
			debug_log.info('qsub prepared 1 ' )
			#DrmaaSessionObj.initialize()
			jt = s.createJobTemplate()
			debug_log.info('qsub prepared 3 ' )
			jt.remoteCommand = FileName
			debug_log.info('qsub prepared 4 ' )
			jt.workingDirectory = os.path.dirname(FileName)
			debug_log.info('qsub prepared 5 ' )
			jt.nativeSpecification = para
			debug_log.info('qsub prepared 6 ' )
		#jt.nativeSpecification = "-b no -shell yes -l vf=1G,p=1,h_vmem=1G -P vip -q sci.q"
			jobid = s.runJob(jt)
			debug_log.info('qsub prepared 7 ' )
			s.deleteJobTemplate(jt)
			debug_log.info('qsub prepared 8 ' )
		debug_log.info('qsub prepared 9 ' )
		return jobid
	except Exception as e :
		print("Exception :{0}".format(e)) 
		debug_log.info("Exception :{0}".format(e))
		return 'error'

class Script():
	def __init__(self , script, count , name):
		self.count = count
		self.name = name 
		self.script = script
		self.queue = ''
		self.up = ''
		self.uvf = ''
		self.uh_vmem = ''
		self.other = ''
		self.max_cycle = ''
		self.jobid = ''
		self.finish = 0
		self.state = 'waiting'
		self.oldp = ''
		self.oldvf=''
		self.oldh_vmem =''
		self.job_obj = ''
		self.qsub_time = 0
		self.qsub_P = ''
	def to_list(self):
		return [self.count , self.name, self.script,  self.queue,     self.up ,    self.uvf,          self.uh_vmem,  self.other,
				self.oldp ,    self.oldvf , self.oldh_vmem, self.max_cycle, self.jobid, self.qsub_time ,  self.finish , self.state]
	def add_attrib(self, queue, p , vf , h_vmem,other, oldp , oldvf , oldhvmem , max_cycle , state , finish , qsub_time ,job_id , qsub_P):
		self.queue = queue
		self.max_cycle = max_cycle
		self.up = p
		self.uvf = vf
		self.uh_vmem = h_vmem
		self.other = other
		self.oldp = oldp
		self.oldvf = oldvf
		self.oldh_vmem = oldhvmem
		self.state = state
		self.finish = finish
		self.qsub_time = qsub_time
		self.jobid = job_id
		self.qsub_P = qsub_P
	def resource(self):
		r_value = ' vf={0.uvf},p={0.up},h_vmem={0.uh_vmem} '.format(self)
		if len(self.other) > 0 :
			for i in self.other : 
				r_value += ',{0}'.format(i)
		return r_value
	def qsub(self):
		para = "-b no -shell yes -l {resource} -q {queue}".format(resource=self.resource(),queue=self.queue)
		if self.qsub_P  and self.qsub_P != 'none':
			para += ' -P {0}'.format(self.qsub_P)
		debug_log.info('qsub ' + self.script + ' ' +para)
		count = 0 
		while count <=5 : 
			count += 1 
			self.jobid = Submit(self.script, para)
			#"check jobid ?"
			if self.jobid == 'error':
				debug_log.info('qsub failed' )
				time.sleep(60)
			else:
				debug_log.info('qsub success ' )
				self.state = 'qsub'
				self.qsub_time += 1
				return self.jobid
		else:
			debug_log.error('Exit: qsub ' + self.script + ' ' +para)
			sys.exit(1)

	def qsub_to_run(self):
		self.time_last = time_current() #time.time()
		self.state = 'running'
	def qhold(self):
		cmd = 'qhold {0.jobid}'.format(self)
		if os.system(cmd):
			pass
		else:
			self.state = 'hold'
			self.job_obj.qhold_time = time_current() #time.time()
	def qrelease(self , db):
		qstat = query(db , 'JOB' , 'JOB_ID' , self.jobid   , '*')
		if len(qstat) > 0 : 
			cmd = 'qrls {0.jobid}'.format(self)
			if os.system(cmd):
				debug_log.warn('qrlease {0} {1} fail'.format(self.name , self.jobid ))
			else:
				self.state = 'running'
				self.job_obj.qrls_time = time_current() #time.time()
		else:
			self.state = 'running'
			self.job_obj.qrls_time = time_current() #time.time()
	def qdel(self , db):
		#cmd = 'qdel {0.jobid}'.format(self)
		self.up = self.up if self.up > int(self.job_obj.uthread + 1 ) else int(self.job_obj.uthread + 1)
		self.uvf = self.job_obj.umaxvmem if float(self.job_obj.umaxvmem) > MEMORY_INCREASE_RATIO * float(self.uvf) else MEMORY_INCREASE_RATIO  * self.uvf
		self.uh_vmem = self.oldh_vmem if self.oldh_vmem > hvmem_ratio * self.uvf else hvmem_ratio  * self.uvf
		r_dict = dict(zip(['OLDP' ,'OLDVF','OLDHVMEM'],[self.oldp,self.oldvf,self.oldh_vmem]))
		db.update('SCRIPT' , r_dict , 'ID', self.count)
		with drmaa.Session() as s:
			try:
				s.control( '{0.jobid}'.format(self) , drmaa.JobControlAction.TERMINATE)
				self.state = 'break'
			except drmaa.errors.InvalidJobException:
				self.state = 'break'
			except Exception as e :
				debug_log.warn("Exception :" + e) 
				debug_log.warn('WARNING : qdel {0.jobid} failed'.format(self))
				sys.exit(e)
		#if os.system(cmd):
		#	debug_log.warn('WARNING : qdel {0.jobid}'.format(self))
		#else:
		#	self.state = 'break'
	def __str__(self):
		return "\t".join( [ str(i) for i in self.to_list()])
	def updatedb(self , db):
		header = ['ID' , 'NAME', 'SCRIPT' , 'QUEUE', 'P','VF','HVMEM' , 'OTHER' , 'OLDP' , 'OLDVF' ,'OLDHVMEM' , 'MAXCYCLE','JOBID','QSUBTIME', 'FINISH' ,'STATE']
		values = self.to_list()
		record = dict(zip(header , values ))
		db.update('SCRIPT', record , 'ID' , self.count)

class Job():
	def __init__(self , jobid):
		self.jobid = jobid
		self.name = ''
		self.script = ''
		self.state = ''
		self.p = ''
		self.vf = ''
		self.h_vmem = ''
		self.other = ''
		self.node = ''
		self.ucpu = 0
		self.umem = 0
		self.uio = 0
		self.umaxvmem = 0  
		self.uthread =  0 
		self.start = ''
		self.end = ''
		self.finish = ''
		self.kill = ''
		self.died = ''
		self.counter_qstat = 0
		self.qhold_time = time_current()
		self.qrls_time = time_current()
		self.time_last = ''
		self.cputime_last = 0 
	def fill_information(self , script_obj):
		self.name = script_obj.name
		self.script = script_obj.script
		self.p = script_obj.up
		self.vf  = script_obj.uvf
		self.h_vmem = script_obj.uh_vmem
		self.other = script_obj.other
	def exceed_resource(self):
		exceed_cpu = self.uthread / self.p
		exceed_mem = self.umaxvmem / self.vf
		#print(exceed_cpu , exceed_mem )
		return ( exceed_cpu , exceed_mem ) 
	def update(self , content):
		self.state = content[2]
		if self.state == 'r':
			self.umaxvmem = self.umaxvmem if self.umaxvmem > float(content[11]) else float(content[11])
			self.ucpu = float(content[9])
			self.umem = float(content[10])
			self.uio  = self.uio if self.uio > float(content[13]) else float(content[13])
			self.uthread = float(content[12])
			self.node = content[4]
	def slow_node(self):
		self.counter_qstat += 1 
		if  self.counter_qstat > 1 and self.counter_qstat % 3 != 0 : 
			return False
		self.time_current = time_current() #time.time()
		#if not self.qstat():return False
		hold_time = (datetime.strptime(self.qrls_time , '%Y-%m-%d %H:%M:%S') - datetime.strptime(self.qhold_time , '%Y-%m-%d %H:%M:%S')).total_seconds()
		if self.state == 'r' :
			if self.time_last == '':
				self.time_last = self.time_current
				#self.cputime_current = self.cputime_last 
			if (datetime.strptime(self.time_current,'%Y-%m-%d %H:%M:%S') - datetime.strptime(self.time_last, '%Y-%m-%d %H:%M:%S')).total_seconds() - hold_time  > 3600 : 
				if self.ucpu - self.cputime_last - hold_time  < 10 :
					return True
				else : 
					self.cputime_last = self.ucpu 
					self.time_last = self.time_current
		else:
			return False
	def to_list(self):
		return [self.jobid , self.name , self.script , self.state ,   self.p ,       self.vf ,    self.h_vmem , self.other , self.node , 
				self.ucpu ,  self.umem , self.uio ,    self.umaxvmem ,self.uthread , self.start , self.end ,    self.finish ,self.kill ,self.died ]
	def __str__(self):
		return "|".join(self.to_list())
	def updatedb(self , db):
		header = ['JOBID' ,'NAME','SCRIPT'  ,'STATE'  ,'P'    ,'VF'   ,'HVMEM' , 'OTHER' , 'NODE'  , 'UCPU' ,
				   'UMEM' ,'UIO' ,'UMAXVMEM','UTHREAD','START','END'  ,'FINISH' ,'KILL'  ,'DIED']
		values = self.to_list()
		record = dict(zip(header , values ))
		db.update('JOB', record , 'JOBID' , self.jobid)
	def insertdb(self , db):
		header = ['JOBID' ,'NAME','SCRIPT'  ,'STATE'  ,'P'    ,'VF'   ,'HVMEM'  , 'OTHER' ,  'NODE'  , 'UCPU' ,
				   'UMEM' ,'UIO' ,'UMAXVMEM','UTHREAD','START','END'  ,'FINISH' ,'KILL'  ,'DIED']
		values = self.to_list()
		db.insert('JOB', header , values)

def popen(cmd , max_count = 10  ):
	count = 0 
	timeout = 60 
	if cmd.find('du') > -1 :
		timeout = 300
	
	while count < max_count :
		child = Command( cmd )
		if child.run( timeout = timeout ) == 0 :
		#child.wait()
		#if child.poll() == 0 :
			return [  i.decode() for i in child.content.split(b'\n') if len(i) ] 
		else:
			#print('{0} is error\n'.format(cmd))
			debug_log.warn('{0} failed , repeat {1} time'.format(cmd , count ))
			#print('{0} failed , repeat {1} time'.format(cmd , count ))
			if cmd.find('du') > -1 :
				return ['0']
			elif cmd.startswith('qstat -j'):
				return 'qstat failed'
			else:
				pass
		count += 1 
		time.sleep(300)
	else :
		debug_log.error('Exit : Command {0} is error'.format(cmd))
		#print('{0} is error\n'.format(cmd))
		if cmd.find('du') > -1 :
			return ['0']
		else:
			sys.exit(1)

def makedir(dir):
	if not os.path.isdir(dir):
		if os.system('mkdir -p {0}'.format(dir)):
			debug_log.error("Exit : cannot mkdir -p {0}".format(dir))
			sys.exit("cannot mkdir -p {0}".format(dir))

def generate_split_shell( args ):
	shell , line_interval ,job_prefix = args.input , args.line ,  args.prefix
	pat = re.compile(';\s*;')
	qsub_dir = '{0}.qsub'.format(os.path.abspath(shell))
	makedir(qsub_dir)
	shell_out = '' 
	job_dict = {}
	job_count = 0 
	with open(shell) as f_in:
		for line_count , content in enumerate(f_in):
			content = content.lstrip().rstrip().replace('&' , ';')
			if line_count % line_interval == 0 :
				if not shell_out == '' :
					shell_out.write('echo This-Work-is-Completed!\n')
					shell_out.close()
					debug_log.info('{0} close'.format(script_file))
				job_count += 1
				job_count_str = '{0:0>4}'.format(job_count)
				script_file = '{0}/{1}_{2}.sh'.format(qsub_dir , job_prefix, job_count_str) 
				shell_out = open(script_file , 'w')
				debug_log.info('{0} open'.format(script_file))
				job_dict[job_count_str] = Script(script_file , job_count_str , '{0}_{1}'.format(job_prefix, job_count_str))
			if not bool(content) : continue
			if content.startswith('#'):continue
			content = content.rstrip(';')
			while(pat.search(content)):
				content = pat.sub(';',content)
			shell_out.write('{0} &&\\\n'.format(content))
		else:
			if shell_out == '' :
				debug_log.error('{0} was empty'.format(shell))
				sys.exit(1)
			else :
				shell_out.write('echo This-Work-is-Completed!\n')
				shell_out.close()
	return job_dict

def is_dir_full(flag_nodu , db_obj , dir , input):
	if flag_nodu: 
		cursor = db_obj.query('PARAMETER' , 'ID'  , '1' ,  ['QUOTA'])
		quota = cursor[0][0]
		#print("dir:{1}  QUOTA:{0} ".format(quota , dir ))
		if quota == '100000000000000000G' :
			return False 
		size =  "".join(popen('du -bs --apparent-size {0}'.format(dir)))
		size = int(size.split()[0])
		debug_log.info("QUOTA:{0} size {1} dir {2}".format(quota , size ,dir))
		qq = transfer(quota)
		if size >0.95*qq:
			try:
				debug_log.info("prepare send directory full email")
				debug_log.info('ssh c0008 2>/dev/null /annoroad/share/software/install/Python-3.3.2/bin/python3 {0}/jump_send_email.py -i {1}.db -full'.format( bindir, os.path.abspath(input)))
				popen(' ssh c0008 /annoroad/share/software/install/Python-3.3.2/bin/python3 {0}/jump_send_email.py -i {1}.db -full -m 0 '.format( bindir, os.path.abspath(input)))
				if os.path.isfile('{0}.db.ini'.format( args.input)):
					popen(' rm {0}.db.ini'.format(input))
				time.sleep(600)
			except:
				pass
			return True
		else:
			return False
	else :
		return False

def query(db , table_name , key , a_value , cols):
	cmd = ''
	if a_value == '*':
		cmd = 'SELECT {3} FROM {0} ; '.format(table_name , key , a_value , cols)
	else:
		cmd = 'SELECT {3} FROM {0} WHERE {1} = "{2}"; '.format(table_name , key , a_value , cols)
	#print(cmd , db)
	conn = sqlite3.connect(db)
	count = 0 
	while ( count < 10):
		try:
			cursor = conn.execute(cmd)
			tt =  cursor.fetchall()
			conn.close()
			return tt
		except sqlite3.OperationalError:
			pass
		count += 1
		time.sleep(30)
	debug_log.error('Exit : cannot access to {0}'.format(db))
	sys.exit(1)

def query_all(db , table_name ):
	cmd = 'SELECT * FROM {0} ; '.format(table_name )
	#print(cmd)
	count = 0 
	conn = sqlite3.connect(db)
	while ( count < 10):
		try:
			cursor = conn.execute(cmd)
			tt =  cursor.fetchall()
			conn.close()
			return tt
		except sqlite3.OperationalError:
			pass
		count += 1
		time.sleep(30)
	debug_log.error('Exit : cannot access to {0}'.format(db))
	sys.exit(1)

def check_die_node(db):
	die_node = query(db , 'NODE' , 'DIED' , 1  , 'NAME')
	return [j for i in die_node for j in i]

def check_o_file(a_script):
	if os.path.isfile('{0.script}.o{0.jobid}'.format(a_script)):
		cmd = 'tail -1 {0.script}.o{0.jobid}'.format(a_script)
		last_line = "".join(popen(cmd))
		time.sleep(3)
		if last_line.rstrip().endswith('This-Work-is-Completed!'):
			return True
		else:
			return False
	else:
		return False

def check_e_file(a_script):
	if os.path.isfile('{0.script}.e{0.jobid}'.format(a_script)):
		cmd = 'tail -1 {0.script}.e{0.jobid}'.format(a_script)
		last_line = "".join(popen(cmd))
		if last_line.rstrip().endswith('MemoryError'):
			return True
		else:
			return False
	else:
		return False


def output_log(logfile , done_job_list , obj_dict ):
	#print("done in" , [i.key for i in done_job_list])
	finish_job_list = []
	with open(logfile , 'a') as f_out:
		total , finish = 0 ,0 
		for i in obj_dict:
			total += 1 
			if obj_dict[i].state == 'done' or obj_dict[i].state == 'finished':
				finish += 1
				finish_job_list.append(obj_dict[i])
		
		if finish > len(done_job_list):
			f_out.write('[Process]:\t{0}/{1} finished\n'.format(finish , total))
			f_out.write('[Finished]:\t{0}\n'.format("\t".join([i.count for i in finish_job_list])))
			done_job_list = finish_job_list
	#print("done out " , [i.key for i in done_job_list])
	return done_job_list

def update_parameters(logfile , db_obj):
	quota , maxjob = '100000000000000000G' , 4 
	with open(logfile) as f_in:
		for i in f_in:
			tmp = i.split()
			if i.startswith('DISK_QUOTA'):
				quota = tmp[1]
			elif i.startswith('Max_Jobs'):
				maxjob = int(tmp[1])
	db_obj.update('PARAMETER',{"QUOTA":quota , "MAXJOBS":maxjob},'ID' , '1')
	return quota , maxjob

def running_count(db , all_script_list) :
	all_running_script = [str(j) for i in  query(db , 'JOB' , 'STATE' , '*', 'JOB_ID') for j in i ]
	query_time_now = query_all(db , 'TIME' )[0][0]
	all_job_id = [i.jobid for i in all_script_list]
	#print(all_job_id)
	count_running_job = set(all_running_script).intersection(all_job_id)
	return(query_time_now , len(count_running_job))

def time_current():
	time_now = datetime.now()
	return time_now.strftime('%Y-%m-%d %H:%M:%S')

def get_unfinish(all):
	r_list = []
	for i in all :
		if i.state in ['done','finished' , 'fail']:
			pass
		else :
			r_list.append(i)
	return r_list

def guard_objs(obj_dict , args ,logfile , db , db_obj):
	all_script_list = [obj_dict[i] for i in sorted(obj_dict)] 
	done_job_list = []
	query_time_last , count_running_job = running_count(db , all_script_list)
	while(all_script_list):
		done_job_list = output_log(logfile , done_job_list , obj_dict)
		dead_node = check_die_node(db)
		if len(dead_node) : 
			debug_log.warn('dead node : {0}'.format(' '.join(dead_node)))
		quota , maxjob = update_parameters(logfile , db_obj)
		query_time_now , count_running_job = running_count(db , all_script_list)
		#print("line579", query_time_now , query_time_last ,  count_running_job )
		debug_log.info('check qstat database , now is  {0} , database time is  {1} , running job count is {2}'.format(query_time_now , query_time_last ,  count_running_job )) 
		if query_time_last == query_time_now :
			time.sleep(30)
			continue
		else:
			query_time_last = query_time_now
		#print('check database new')
		bool_dir_full = is_dir_full( args.nodu , db_obj , args.analysis_dir ,args.input )
		#print('direcotry full {0}'.format(bool_dir_full))
		debug_log.info('direcotry full {0}'.format(bool_dir_full))
		if bool_dir_full:
			for a_script in all_script_list:
				if a_script.state in ['qsub' ,'running']:
					a_script.qhold()
					a_script.job_obj.qhold_time = time_current()
					a_script.job_obj.updatedb(db_obj)
					a_script.updatedb(db_obj)
			time.sleep(args.interval)
			continue 
		else : 
			for a_script in all_script_list:
				if (count_running_job >=  maxjob): 
					break
				bool_finish =  check_o_file(a_script)
				bool_memory_error = check_e_file(a_script)
				#print(a_script.name , a_script.state , "check o",  bool_finish)
				debug_log.info('{0} {1} check .o file {2}'.format(a_script.name , a_script.state ,  bool_finish))
				if a_script.state in ['hold']:
					a_script.qrelease(db)
					a_script.job_obj.qrls_time = time_current()
					a_script.job_obj.updatedb(db_obj)
					a_script.updatedb(db_obj)
				elif a_script.state in ['waiting' , 'break']:
					if bool_finish :
						a_script.state = 'finished'
						a_script.finish = 1
						a_script.job_obj.finish = 1
						a_script.job_obj.end = time_current()
						a_script.job_obj.updatedb(db_obj)
						a_script.updatedb(db_obj)
					else:
						if a_script.qsub_time >= a_script.max_cycle : 
							a_script.state = 'fail'
							a_script.job_obj.updatedb(db_obj)
							a_script.updatedb(db_obj)
						else:
							job_id = a_script.qsub()     #remove logfile 20180209
							#print(a_script.name , 'qsub success')
							debug_log.info('{0} qsub success'.format(a_script.name))
							a_script.job_obj = Job(job_id)
							a_script.job_obj.start = time_current()
							a_script.job_obj.fill_information(a_script)
							a_script.job_obj.insertdb(db_obj)
							a_script.updatedb(db_obj)
							count_running_job += 1 
				elif a_script.state in ['qsub' , 'running']:
					qstat = query(db , 'JOB' , 'JOB_ID' , a_script.jobid   , '*')
					if len(qstat) > 0 :
						#print(qstat[0])
						a_job = a_script.job_obj
						a_job.update( qstat[0] )
						if a_job.state == 'r': 
							bool_kill_flag =  True
							if a_job.node in dead_node:
								debug_log.warn('WARNING : job {0} on dead node {1}'.format(a_job.jobid , dead_node))
							#elif a_job.slow_node() :
							#	debug_log.warn( 'WARNING :job {0} on slow node {1}'.format( a_job.jobid , a_job.slow_node() ))
							elif a_job.exceed_resource()[0] > args.kill_cpu:
								debug_log.warn( 'WARNING:job {0} exceed cpu {1}'.format( a_job.jobid , a_job.exceed_resource()[0]  ))
							elif a_job.exceed_resource()[1] >args.kill_mem :
								debug_log.warn(' WARNING:job {0} exceed mem {1}'.format( a_job.jobid, a_job.exceed_resource()[1] ))
							else :
								bool_kill_flag =  False

							if bool_kill_flag : 
								a_script.qdel(db_obj)
								a_job.end = time_current()
								a_job.kill = 1 
								count_running_job -= 1
								a_script.state = 'break'
								a_script.jobid = ''
								a_script.job_obj.updatedb(db_obj)
								a_script.updatedb(db_obj)
							else:
								if a_script.state == 'qsub':
									a_script.qsub_to_run()
									a_job.start = time_current()
									a_script.job_obj.updatedb(db_obj)
									a_script.updatedb(db_obj)
								else : 
									a_script.job_obj.updatedb(db_obj)
									a_script.updatedb(db_obj)
						elif a_job.state in ['hqw' , 'hr' , 'ht' ]: 
							a_script.state = 'hold'
							a_script.job_obj.updatedb(db_obj)
							a_script.updatedb(db_obj)
						elif a_job.state in ['qw', 't']:
							if a_script.state == 'qsub':
								pass
							else : 
								#print(a_script.name , 'error')
								debug_log.info('{0} error'.format(a_script.name))
						else :
							#print(a_job.__str__())
							debug_log.warn('WARNING : {0} was qdel'.format(a_job.jobid))
							a_script.qdel(db_obj)
							a_job.end = time_current()
							a_job.kill = 1
							count_running_job -= 1
							a_script.state = 'break'
							a_script.job_obj.updatedb(db_obj)
							a_script.updatedb(db_obj)
							a_script.jobid = ''
					else:
						if bool_finish:
							a_script.state = 'finished'
							a_script.finish = 1
							a_script.job_obj.end = time_current()
							a_script.job_obj.finish = 1
							a_script.job_obj.updatedb(db_obj)
							a_script.updatedb(db_obj)
						else :
							if a_script.qsub_time >= a_script.max_cycle : 
								a_script.state = 'fail'
							else : 
								a_script.state = 'break'
							if bool_memory_error:
								a_script.up = a_script.up if a_script.up > int(a_script.job_obj.uthread + 1) else int(a_script.job_obj.uthread + 1 )
								a_script.uvf = a_script.job_obj.umaxvmem if float(a_script.job_obj.umaxvmem) > MEMORY_INCREASE_RATIO * float(a_script.uvf) else  MEMORY_INCREASE_RATIO *float(a_script.uvf)
								a_script.uh_vmem = a_script.oldh_vmem if a_script.oldh_vmem > hvmem_ratio * a_script.uvf else hvmem_ratio * a_script.uvf
								#a_script.uvf = a_script.oldh_vmem
								#a_script.uh_vmem = 1.5 * a_script.oldh_vmem
							else :
								a_script.up = a_script.up if a_script.up > int(a_script.job_obj.uthread + 1) else int(a_script.job_obj.uthread + 1 )
								a_script.uvf = a_script.job_obj.umaxvmem if float(a_script.job_obj.umaxvmem) > MEMORY_INCREASE_RATIO * float(a_script.uvf) else MEMORY_INCREASE_RATIO * float(a_script.uvf)
								a_script.uh_vmem = a_script.oldh_vmem if a_script.oldh_vmem > hvmem_ratio * a_script.uvf else hvmem_ratio * a_script.uvf
								#self.uvf = self.job_obj.umaxvmem if float(self.job_obj.umaxvmem) > float(self.oldvf) else  self.oldvf
								#self.uh_vmem = self.oldh_vmem if self.oldh_vmem > 1.5 * self.uvf else 1.5 * self.uvf
							a_script.job_obj.end = time_current()
							a_script.job_obj.updatedb(db_obj)
							a_script.updatedb(db_obj)
				else : 
					pass
		#print("sleep start")
		debug_log.info("sleep start")
		time.sleep(args.interval)
		all_script_list = get_unfinish(all_script_list)
		#print("sleep finish")
		debug_log.info("sleep finish")
		if args.noreqsub : break 
	else :
		done_job_list = output_log(logfile , done_job_list , obj_dict)

def output_finish_log(obj_dict , logfile):
	all_finish_flag = True
	with open(logfile ,'a') as f_out:
		failed_list = []
		for i in obj_dict:
			a_script = obj_dict[i]
			if a_script.state == 'finished':
				f_out.write(a_script.__str__() + '\n')
			elif a_script.state == 'fail':
				failed_list.append(a_script.name)
			else:
				debug_log.info("finish {0.name} {0.state}".format(a_script))
		if len(failed_list) > 0:
			f_out.write('{0} is not finish\n'.format("\t".join(failed_list)))
			debug_log.error('Exit : {0} is not finish'.format("\t".join(failed_list)))
			all_finish_flag = False
		else:
			f_out.write('All jobs finished!')
			debug_log.info('All jobs finished!')
			#sys.exit(0)
	return all_finish_flag

def parse_resource(resource):
	p=1
	vf = 1024**3
	h_vmem_flag  = False 
	h_vmem = 0
	other = []
	for mm in re.split('-l|,' , resource):
		mm = mm.strip()
		if mm.find('=') > 0 :
			tt = mm.split('=')
			if tt[0] == 'vf':
				vf = transfer(tt[1])
			elif tt[0] == 'p' :
				p = float(tt[1])
			elif tt[0] == 'h_vmem':
				h_vmem = transfer(tt[1])
				h_vmem_flag = True 
			else :
				other.append(mm)
	if not h_vmem_flag : 
		h_vmem = hvmem_ratio * vf
	return p , vf ,h_vmem , ",".join(other)

def transfer(count):
	c_maxvmem = 0
	if count.endswith('T'):
		c_maxvmem = float(count.replace('T', '')) * 1024 * 1024 * 1024 * 1024
	elif count.endswith('t'):
		c_maxvmem = float(count.replace('t', '')) * 1024 * 1024 * 1024 * 1024
	elif count.endswith('G'):
		c_maxvmem = float(count.replace('G', '')) * 1024 * 1024 * 1024
	elif count.endswith('g'):
		c_maxvmem = float(count.replace('g', '')) * 1024 * 1024 * 1024
	elif count.endswith('M'):
		c_maxvmem = float(count.replace('M', '')) * 1024 * 1024 
	elif count.endswith('m'):
		c_maxvmem = float(count.replace('m', '')) * 1024 * 1024
	elif count.endswith('K'):
		c_maxvmem = float(count.replace('K', '')) * 1024
	elif count.endswith('k'):
		c_maxvmem = float(count.replace('k', '')) * 1024
	elif count == 'N/A':
		c_maxvmem = 0
	else:
		try:
			c_maxvmem = float(count)
		except:
			print("{0} is error at 894".format(count))
			c_maxvmem = 1024 * 1024 * 1024
	return c_maxvmem        
	

def check_previous_job(f_in):
	flag_pre = False
	for a_file in glob.glob('{0}.*.log'.format(f_in)):
		pid = a_file.split(r'.')[-2]
		if pid.isnumeric() : 
			tt = int(pid)
			try:
				os.kill(tt , 0)
			except OSError:
				pass
			else:
				print("{0} is not finish , please kill it before resubmit".format(tt))
				debug_log.error("Error:{0} is not finish , please kill it before resubmit".format(tt))
				flag_pre = True
	if flag_pre :
		sys.exit()


def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument(help='input file',dest='input')
	#parser.add_argument('-o','--output',help='output log  file',dest='output',type=argparse.FileType('w'),required=True)
	parser.add_argument('-l','--lines',  help='line number , default is [1] ',dest='line', type=int , default=1)
	parser.add_argument('-d','--database',  help='database , default is /annoroad/data1/bioinfo/PMO/liutao/database/sge.db ',dest='database', default="/annoroad/data1/bioinfo/PMO/liutao/database/sge.db")
	parser.add_argument('-m','--maxjob',help='max job number , default is [4]', dest='maxjob', type=int ,default=4)
	parser.add_argument('-i','--interval',help='interval check time , default is [180]',dest='interval', type=int , default=180)
	parser.add_argument('-q','--queue',help='job queue , default is [sci.q]' , dest='queue', default='sci.q')
	parser.add_argument('-nr','--noreqsub',help='do not reqsub failed job ,default is reqsub', dest='noreqsub', action='store_true')
	parser.add_argument('-nc','--nocontinue',help='do not continue with unfinish log , default is continue',dest='nocontinue',action='store_true')
	parser.add_argument('-re','--resource',help='resouce list ,default is [ "vf=1G,p=1,h_vmem=1G" ] ', dest = 'resource',default=' vf=1G,p=1,h_vmem=1G ')
	parser.add_argument('-P','--P',help='project in qsub  ', dest = 'qsub_P',default='none')
	parser.add_argument('-kc','--killcpu',help='kill job when thread exceed',dest='kill_cpu' , default= 1.3)
	parser.add_argument('-km','--killmem',help='kill job when mem exceed',dest='kill_mem' ,  default= 1 )
	parser.add_argument('-hv','--hvmem',help='hvmem',dest='hvmem',default= 1)
	parser.add_argument('-prefix','--jobprefix', help='job prefix ,default is [work]',dest='prefix' ,default='work')
	parser.add_argument('-maxcycle','--maxcycle',help='max cycle , ,default is [3]', dest='max_cycle',default=3,type=int)
	parser.add_argument('-quota', '--quota', help='quota,default is [100000000000000000G]', dest='quota' , default='100000000000000000G')
	parser.add_argument('-analysis_dir', '--analysis_dir', help='analysis dir,default is [shell/../..]',dest='analysis_dir' )
	parser.add_argument('-nodu','--nodu',help='do not check disk,default is [du]',dest='nodu',action='store_false')
	parser.add_argument('-js','--jobidstart',help='job id start number',dest='pre_job_count',default= 0 ,type = int)
	parser.add_argument('-mail','--mail_retry' , help='max retry times',dest='retry', default = 0 , type=int)
	
	#sys.exit()
	home_dir =  os.environ['HOME']
	if not os.path.isfile('{0}/.email/.email.txt'.format(home_dir)):
		print('''please config your email address first in ~/.email/.email.txt:
[DEFAULT]
Max_count = 5 
Sleep_time = 60 
[HEADER]
Addressor = **@genome.cn
Password = **** 
Receiver = chengfangtu@annoroad.com; 
Copy = yanxunsu@annoroad.com 
Server = smtp.exmail.qq.com 
Receive_server = pop.exmail.qq.com
''')
		sys.exit(1)
	args=parser.parse_args()

	global hvmem_ratio
	hvmem_ratio = args.hvmem

	work_dir = os.path.abspath(os.path.dirname(args.input))
	if args.analysis_dir == None:
		args.analysis_dir = '{0}/../../'.format(os.path.abspath(work_dir))
	print(args.noreqsub , args.nocontinue , args.nodu , args.analysis_dir)
	global debug_log
	debug_log = mylogger(args.input)
	debug_log.info('{0}.db'.format(args.input))
	if args.nocontinue:
		debug_log.info('rm -rf {0}.db {0}.qsub'.format(args.input))
		popen('rm -rf {0}.db {0}.qsub'.format(args.input))
	db_obj = Database( '{0}.db'.format(args.input))
	bool_exists = db_obj.check_table_exists()
	
	check_previous_job(args.input)
	logfile = '{0}.{1}.log'.format(args.input , os.getpid())


	with open(logfile,'w') as f_out:
		f_out.write('DISK_QUOTA\t{0}\nMax_Jobs\t{1}\n'.format( args.quota, args.maxjob))
	
	parameter_table_header = ['ID' ,'MAXJOBS' , 'INPUT', 'QUOTA']
	script_table_header = ['ID' , 'NAME', 'SCRIPT' , 'QUEUE', 'P','VF','HVMEM' , 'OTHER' , 'OLDP' , 'OLDVF' ,'OLDHVMEM' , 'MAXCYCLE','JOBID','QSUBTIME', 'FINISH' ,'STATE']
	job_table_header = ['JOBID' ,'NAME','SCRIPT'  ,'STATE'  ,'P'    ,'VF'   ,'HVMEM' ,'OTHER' ,'NODE'  , 'UCPU' ,
				        'UMEM' ,'UIO' ,'UMAXVMEM','UTHREAD','START','END'  ,'FINISH' ,'KILL'  ,'DIED']
	p , vf , h_vmem , other = parse_resource(args.resource)
	

	obj_dict = {}
	if not bool_exists:
		db_obj.delete_all_content( 'PARAMETER' )
		db_obj.insert( 'PARAMETER' , parameter_table_header , ['1',args.maxjob, args.input, args.quota])
		
		obj_dict  = generate_split_shell( args )
		for job_count_str in sorted(obj_dict) :
			obj_dict[job_count_str].add_attrib(args.queue ,p,vf,h_vmem,other, p , vf , h_vmem , args.max_cycle ,'waiting' , 0 , 0  , '' ,args.qsub_P)
			db_obj.insert('SCRIPT' , script_table_header , obj_dict[job_count_str].to_list())
	else :
		for i in db_obj.query_all('SCRIPT'):
			(ID,NAME,SCRIPT,QUEUE,P,VF,HVMEM,OTHER,OLDP,OLDVF,OLDHVMEM,MAXCYCLE,JOBID,QSUBTIME,FINISH,STATE) = i
			a_script = Script(SCRIPT , ID , NAME)
			obj_dict[ID] = a_script
			job_id = JOBID
			if (p , vf , h_vmem , other) == (float(P) , float(VF) , float(HVMEM) , OTHER ):
				a_script.add_attrib(QUEUE , float(P) , float(VF) , float(HVMEM) , OTHER , float(OLDP) , float(OLDVF) , float(OLDHVMEM) ,  args.max_cycle ,STATE , int(FINISH) , 0  , job_id , args.qsub_P)
			else:
				a_script.add_attrib(QUEUE , p , vf , h_vmem , other , p , vf , h_vmem , args.max_cycle ,'waiting' , 0, 0  , '' , args.qsub_P)
			if job_id:
				a_script.jobid = job_id
				a_script.job_obj = Job(job_id)
				a_script.job_obj.fill_information(a_script)
				qstat = query(args.database , 'JOB' , 'JOB_ID' , job_id   , '*')
				bool_finish =  check_o_file(a_script)
				#if a_script.state in []:
				#print(a_script.state , a_script.name , 'check o', bool_finish)
				debug_log.info('{0} {1} check .o file {2}'.format(a_script.name , a_script.state ,  bool_finish))
				if len(qstat) >  0 :
					a_script.job_obj.updatedb(db_obj)
					a_script.updatedb( db_obj )
				else:
					if bool_finish:
						a_script.state = 'finished'
						a_script.finish = 1
						a_script.job_obj.finish = 1
						a_script.job_obj.updatedb(db_obj)
						a_script.updatedb( db_obj )
					else :
						a_script.job_obj = ''
						a_script.state =  'waiting'
						a_script.updatedb( db_obj )
	debug_log.info("guard job start")
	guard_objs(obj_dict , args , logfile , args.database ,db_obj)
	all_finish_flag = output_finish_log(obj_dict , logfile)
	db_obj.close()

	try:
		debug_log.info("prepare send email")
		debug_log.info(' ssh c0008 2>/dev/null /annoroad/share/software/install/Python-3.3.2/bin/python3 {0}/jump_send_email.py -i {1}.db -m  {2}'.format( bindir, os.path.abspath(args.input) , args.retry ))
		return_info = popen(' ssh c0008 /annoroad/share/software/install/Python-3.3.2/bin/python3 {0}/jump_send_email.py -i {1}.db -m {2} '.format( bindir, os.path.abspath(args.input , args.retry)))
		debug_log.info("|".join(return_info))
		debug_log.info('send email success')
	except:
		pass
	if os.path.isfile('{0}.db.ini'.format( args.input)):
		popen(' rm {0}.db.ini'.format( args.input))
	if all_finish_flag :
		sys.exit(0)
	else :
		sys.exit(1)

if __name__ == '__main__':
	main()
