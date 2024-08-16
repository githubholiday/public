#! /usr/bin/env python3
import argparse
import sys
import os
import re
import xml.etree.ElementTree as xml
from subprocess import check_output
import sqlite3
from datetime import datetime
bindir = os.path.abspath(os.path.dirname(__file__))


__author__='Liu Tao'
__mail__= 'taoliu@annoroad.com'

pat1=re.compile('^\s*$')

class Node():
	def __init__(self , name ) :
		self.name = name 
	def set(self , tag , value ):
		if tag in self.__dict__: print('error')
		self.__setattr__(tag , value )
	def is_died(self):
		self.died = False
		if "down" in self.state : self.died = True
		if "drain" in self.state : self.died = True
		if "drng" in self.state : self.died = True
		if "fail" in self.state : self.died = True
		return self.died
	def to_sql(self):
		pass
	def __str__(self):
		pass


class Job():
	def __init__(self):
		pass
	def set(self , tag , value ):
		if tag in self.__dict__: print('tag error')
		self.__setattr__(tag , value )
	def modify(self , tag , value ):
		self.__setattr__(tag , value )
	def to_sql(self):
		pass
	def __str__(self):
		pass
	def to_list(self , time_now):
		#project_table_header = ['JOB_ID' ,'PERSON' ,'STATE', 'QUEUE','NODE','NUM_PROC','VF','START','CPU','VMEM', 'MAXVMEM','THREAD', 'IO' , 'VEM']
		print(self.__dict__)
		if self.STATE in ['r']:
			self.queue = self.QUEUE
			self.node = self.NODE
			second = (time_now - datetime.strptime( self.START , '%Y-%m-%dT%H:%M:%S')).total_seconds()
			self.cpu = getattr( self , 'NUM_PROC' , 0 )  
			self.vmem = getattr( self , 'VMEM' , 0 )  
			self.maxvmem = getattr( self , 'MAXVMEM' , 0 )  
			self.io = getattr( self , 'IO' , 0 )  
			self.thread = getattr( self , 'THREAD' , 0 )
			#self.thread = float(self.cpu) / second
		else:
			self.queue = '-'
			self.node = '-'
			self.JAT_start_time = '-'
			self.cpu = '-'
			self.vmem = '-'
			self.maxvmem = '-'
			self.thread = '-'
			self.io = '-'

		self.num_proc = getattr( self , 'num_proc' , 1 )  
		self.h_vmem = getattr( self , 'VF' , '1073741824' )  
		self.virtual_free = getattr( self , 'VF' , '1073741824')

		return [self.JOB_ID , self.PERSON , self.STATE , self.queue , self.node , self.num_proc , self.virtual_free , 
				self.h_vmem , self.START , self.cpu , self.vmem , self.maxvmem , self.thread , self.io , self.SCRIPT]


class Database():
	def __init__(self, path ):
		self.database = path
		self.conn = sqlite3.connect( self.database )
	def check_table_exists(self):
		cmd = 'select name from sqlite_master where type = "table";'
		cursor = self.conn.execute(cmd)
		tables = [ j for i in cursor.fetchall() for j in i]
		r_value = True
		print(tables)
		for i in [ 'NODE' , 'JOB' , 'TIME' ]:
			if len(tables) ==  0 or  not i in tables:
				self.create_database( i )
				r_value = False
		return r_value
	def create_database(self , type ):
		cmd = ''
		if type == 'TIME':
			cmd = '''CREATE TABLE TIME
				(NOW     TEXT      PRIMARY KEY ) ; '''
		elif type == 'NODE':
			cmd = '''CREATE TABLE NODE
				(NAME     TEXT      PRIMARY KEY  ,
				 DIED     INT       NOT NULL) ; '''
		elif type == 'JOB':
			cmd = '''CREATE TABLE JOB 
				(JOB_ID   INT       PRIMARY KEY  , 
				 PERSON   TEXT      NOT NULL , 
				 STATE    TEXT      NOT NULL , 
				 QUEUE    TEXT      NOT NULL , 
				 NODE     TEXT      NOT NULL , 
				 NUM_PROC INT       NOT NULL ,
				 VF       TEXT      NOT NULL , 
				 HVMEM       TEXT      NOT NULL , 
				 START    TEXT   ,
				 CPU      TEXT   ,
				 VMEM     TEXT   ,
				 MAXVMEM  TEXT   ,
				 THREAD  TEXT   ,
				 IO      TEXT   ,
				 SCRIPT   TEXT  ); '''
		self.conn.execute(cmd)
	def query(self , table_name , key , a_value , cols):
		cmd = 'SELECT {3} FROM {0} WHERE {1} = "{2}"; '.format(table_name , key , a_value , ",".join(cols))
		print(cmd)
		self.conn.execute(cmd)
	def insert(self , table_name , name_list , value_list):
		names = ",".join(name_list) 
		values = ",".join(['"' + str(i) + '"' for i in  value_list])
		cmd = 'INSERT INTO {0} ({1}) VALUES ({2});'.format(table_name , names , values)
		print(cmd)
		self.conn.execute(cmd)
	def delete(self , table_name , key , a_value):
		cmd = 'DELETE FROM {0} WHERE {1} = {2} ; '.format(table_name , key , a_value)
		print(cmd)
		self.conn.execute(cmd)
	def delete_all_content(self , table_name ):
		cmd = 'DELETE  FROM {0} ;  '.format(table_name )
		print(cmd)
		self.conn.execute(cmd)
	def update(self , table_name , record , key , a_value):
		cmd = ''
		for i , j in record.items():
			cmd += 'UPDATE {0} set {1} = "{2}" WHERE {3} = "{4}";\n'.format(table_name , i, j , key , a_value )
		print(cmd)
		self.conn.execute(cmd)
	def commit(self):
		self.conn.commit()
	def close(self):
		self.conn.close()

def my_popen(command):
	out = check_output(command , shell = True)
	return out 

def qhost():
	command = "sinfo -N"
	result =  my_popen( command )
	node_list = []

	for line in result.decode().split('\n'):
		if line.startswith('NODELIST'):
			continue
		elif line.startswith('-----') or pat1.match(line):
			continue
		else:
			tmp = line.split()
			name , state = tmp[0] , tmp[-1]
			obj_node = Node(name)
			obj_node.set('state' , state)
			node_list.append(obj_node)
	return node_list

def qstat():
	person_list_file = "{0}/person.list".format(bindir)
	if os.path.exists(person_list_file):
		all_person = []
		with open("{0}/person.list".format(bindir)) as fin:
			for line in fin:
				all_person.append(line.rstrip())

		#command = '/opt/gridengine/bin/lx26-amd64/qstat -u \* -xml'
		command = 'squeue -o "%18i %9P %12j %12u %12T %12M %16l %6D %R " -u {0}'.format(",".join(all_person))
	else:	
		#command = '/opt/gridengine/bin/lx26-amd64/qstat -u \* -xml'
		user = os.environ['USER']
		command = 'squeue -o "%18i %9P %12j %12u %12T %12M %16l %6D %R " -u {0}'.format(user)

	print(command)

	result =  my_popen( command )
	'''
	JOBID              PARTITION NAME         USER         STATE        TIME         TIME_LIMIT       NODES  NODELIST(REASON)
	7811846            xahcnorma test         acdgo9idhi   RUNNING      1:14         333-08:00:00     1      c01r2n20
	'''
	job_dict  = {}

	for line in result.decode().split('\n'):
		line = line.strip()
		print(line)
		if line.startswith('JOBID'):
			continue
		elif line.startswith('-----') or pat1.match(line):
			continue
		else :
			tmp = line.split()
			job_id , partition , name , user , state , time , time_limit , nodes , node_list = tmp[0] , tmp[1] , tmp[2] , tmp[3] , tmp[4] , tmp[5] , tmp[6] , tmp[7] , tmp[8]

			if state == 'RUNNING':
				r_state = 'r'
			elif state == 'PENDING':
				r_state = 'qw'
			elif state == 'COMPLETED':
				r_state = 'd'
			elif state == 'CANCELLED':
				r_state = 'c'
			elif state == 'FAILED':
				r_state = 'f'
			elif state == 'TIMEOUT':
				r_state = 't'
			elif state == 'OUT_OF_MEMORY':
				r_state = 'oom'

			obj_job = Job()
			obj_job.set('JOB_ID' , job_id)
			obj_job.set('PERSON' , user)
			obj_job.set('STATE' , r_state)
			obj_job.set('QUEUE' , partition)
			obj_job.set('NODE' , node_list)
		job_dict[ job_id ] = obj_job



	#sacct --format="JobID%20,JobName,User,Partition,NodeList,Elapsed,TotalCPU,State,AllocTRES%50,CPUTime,MaxRSS" -j 7811846"
	#https://rc.byu.edu/wiki/?id=Using+sacct#:~:text=You%20can%20increase%20or%20decrease%20the%20number%20of,name%20column%20to%20have%2040%20characters%2C%20left%20justified%3A

	if len(job_dict) == 0:
		return job_dict
		
	command = 'sacct --format="JobID%20,JobName,User,Partition,NodeList,Elapsed,TotalCPU,State,AllocTRES%50,CPUTime,MaxRSS,AveRSS,ReqMem,MaxDiskRead,MaxDiskWrite,ExitCode%10,NCPUs,Start,Command%10000" -p -j ' + ",".join( job_dict.keys())

	print(command)
	result =  my_popen( command )
	#print(result.decode() , file = sys.stdout)

	for line in result.decode().split('\n'):
		line = line.strip()
		if line.startswith('JobID'):
			continue
		elif line.startswith('-----') or pat1.match(line):
			continue
		else :
			tmp = line.split("|")
			print(tmp)
			job_id , job_name , user , partition , node_list , elapsed , total_cpu , state , alloc_tres , cpu_time , max_rss , ave_rss , req_mem , max_disk_read , max_disk_write , exit_code , ncpus , start , command = tmp[0] , tmp[1] , tmp[2] , tmp[3] , tmp[4] , tmp[5] , tmp[6] , tmp[7] , tmp[8] , tmp[9] , tmp[10] , tmp[11] , tmp[12] , tmp[13] , tmp[14] , tmp[15] , tmp[16] , tmp[17] , " ".join(tmp[18:])
			
			request_cpu=0 
			for a_tres in alloc_tres.split(','):
				if 'cpu' in a_tres:
					request_cpu = int(a_tres.split('=')[1])

			if req_mem.endswith('c'):
				req_mem = req_mem[:-1]
			if req_mem.endswith('n'):
				req_mem = req_mem[:-1]
			
			if not ave_rss: ave_rss = 0
			if not max_rss:	max_rss = 0
			if not max_disk_read: max_disk_read = 0
			if not max_disk_write: max_disk_write = 0
			if not command: command = '-'

			if job_id in job_dict:
				obj_job = job_dict[ job_id ]
				obj_job.set('NUM_PROC' , request_cpu)
				obj_job.set('VF' , req_mem)
				obj_job.set('HVMEM' , req_mem)
				obj_job.set('START' , start)

				obj_job.set('CPU' , ncpus)
				obj_job.set('VMEM' , ave_rss)
				obj_job.set('MAXVMEM' , max_rss)
				obj_job.set('THREAD' , ncpus)
				obj_job.set('IO' , max_disk_read + max_disk_write)
				obj_job.set('SCRIPT' , command)
				job_dict[ job_id ] = obj_job
			else:
				print(job_id , 'not in job_dict')
	return job_dict

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-d','--database', dest='database', 
					 help='database , default is {0}/database/slurm.db'.format(os.environ['HOME']), 
					 default='{0}/database/slurm.db'.format(os.environ['HOME']))
	#parser.add_argument('-m','--mm',help='output file',dest='mm',action='store_false')
	args=parser.parse_args()
	
	db = Database( args.database )
	qhost_result = qhost()
	qstat_result = qstat()
	time_now = datetime.now()
	project_table_header = ['JOB_ID' ,'PERSON' ,'STATE', 'QUEUE','NODE','NUM_PROC','VF','HVMEM' , 'START','CPU','VMEM', 'MAXVMEM','THREAD' , 'IO' , 'SCRIPT']
	if db.check_table_exists():
		db.delete_all_content('TIME')
		db.delete_all_content('NODE')
		db.delete_all_content('JOB')
	db.insert('TIME' , ['NOW'] , [ time_now.strftime('%Y-%m-%d %H:%M:%S')])
	for a_node in qhost_result : 
		db.insert('NODE' , ['NAME' , 'DIED'] , [a_node.name , int(a_node.is_died())]) 
	for a_job in qstat_result : 
		db.insert('JOB' , project_table_header , qstat_result[a_job].to_list( time_now ))
	db.commit()
	db.close()
	print("all job finish")


if __name__ == '__main__':
	main()
