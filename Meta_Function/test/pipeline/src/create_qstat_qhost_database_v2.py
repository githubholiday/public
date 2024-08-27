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

pat1=re.compile('^\s+$')

class Node():
	def __init__(self , name ) :
		self.name = name 
	def set(self , tag , value ):
		if tag in self.__dict__: print('error')
		self.__setattr__(tag , value )
	def is_died(self):
		self.died = False
		for i in ['load_avg' , 'mem_total' , 'mem_used' , 'swap_total' , 'swap_used']:
			if self.__getattribute__(i).find(r'-') > -1 :
				self.died =True
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
		if self.state in ['r']:
			self.queue = self.queue_name.split('@')[0]
			self.node = self.queue_name.split('@')[1]
			second = (time_now - datetime.strptime( self.JAT_start_time , '%Y-%m-%dT%H:%M:%S')).total_seconds()
			self.cpu = getattr( self , 'cpu' , 0 )  
			self.vmem = getattr( self , 'vmem' , 0 )  
			self.maxvmem = getattr( self , 'maxvmem' , 0 )  
			self.io = getattr( self , 'io' , 0 )  
			self.thread = float(self.cpu) / second
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
		self.h_vmem = getattr( self , 'h_vmem' , '1073741824' )  
		self.virtual_free = getattr( self , 'virtual_free' , '1073741824')

		return [self.JB_job_number , self.JB_owner , self.state , self.queue , self.node , self.num_proc , self.virtual_free , 
				self.h_vmem , self.JAT_start_time , self.cpu , self.vmem , self.maxvmem , self.thread , self.io , self.script]


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
	command = "/opt/gridengine/bin/lx26-amd64/qhost -xml"
	result =  my_popen( command )
	root = xml.fromstring( result )
	node_list = []
	for child in root.findall('host'):
		name = child.get('name')
		obj_node = Node(name)
		for i in ['load_avg' , 'mem_total' , 'mem_used' , 'swap_total' , 'swap_used']:
			obj_node.set(i , child.find('.//hostvalue[@name="{0}"]'.format(i)).text)
		node_list.append(obj_node)
	return node_list

def qstat():
	all_person = []
	with open("{0}/person.list".format(bindir)) as fin:
		for line in fin:
			all_person.append(line.rstrip())

	#command = '/opt/gridengine/bin/lx26-amd64/qstat -u \* -xml'
	command = '/opt/gridengine/bin/lx26-amd64/qstat -u {0} -xml'.format(" ".join(all_person))
	result =  my_popen( command )
	root = xml.fromstring( result.decode().replace('\b' , '') )
	job_dict  = {}
	for child in root.findall('.//job_list'):
		job_id = child.find('.//JB_job_number').text
		obj_job = Job()
		for next_child in child.iter():
			obj_job.set( next_child.tag , next_child.text)
		job_dict[ job_id ] = obj_job
	
	command = '/opt/gridengine/bin/lx26-amd64/qstat -xml -j ' + ",".join( job_dict.keys())
	result =  my_popen( command )
	print(result.decode() , file = sys.stdout)
	root = xml.fromstring( result.decode().replace('<>' , '<a>').replace('</>' , '</a>') )
	for child in root.findall('.//djob_info/element'):
		job_id = child.find('.//JB_job_number').text
		if child.find('.//JB_script_file') == None : 
			print(job_id)
			#raise TypeError
			#sys.exit()
			job_script_file = '-' 
		else :
			job_script_file = child.find('.//JB_script_file').text
		if child.find('.//JB_cwd') == None : 
			jb_cwd = ''
		else :
			jb_cwd = child.find('.//JB_cwd').text
		script = job_script_file
		if not job_script_file.startswith('/'):
			script = jb_cwd + '/' + job_script_file
		a_job = job_dict[ job_id ]
		a_job.set('script' , script)
		for next_child in child.findall('.//JB_hard_resource_list/qstat_l_requests'):
			name = next_child.find('.//CE_name').text
			value = next_child.find('.//CE_doubleval').text
			a_job.set(name , value)
		for next_child in child.findall('.//JAT_scaled_usage_list/scaled'):
			name = next_child.find('.//UA_name').text
			value = next_child.find('.//UA_value').text
			a_job.set(name , value)
	for child in root.findall('.//a/ST_name'):
		print(child.tag , child.attrib , child.text)
		job_id = child.text
		a_job = job_dict[ job_id ]
		a_job.modify('state' ,'finish')
	return job_dict

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-d','--database',help='database file',dest='database',default='/annoroad/data1/bioinfo/PMO/liutao/database/sge.db')
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
