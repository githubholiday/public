#! /usr/bin/env python3
# -*- coding: utf-8 -*-  
import argparse
import sys
import os
import re
import sqlite3
bindir = os.path.abspath(os.path.dirname(__file__))
import myemail

__author__='Liu Tao'
__mail__= 'taoliu@annoroad.com'

pat1=re.compile('^\s+$')

def query_all(db , table_name):
	cmd = 'SELECT * FROM {0} ; '.format( table_name )
	conn = sqlite3.connect(db)
	cursor = conn.execute(cmd)
	tt =  cursor.fetchall()
	conn.close()
	return tt

def query(db , table_name , key , value):
	cmd = 'SELECT * FROM {0} WHERE {1} = "{2}" ; '.format( table_name , key , value )
	conn = sqlite3.connect(db)
	cursor = conn.execute(cmd)
	tt =  cursor.fetchall()
	conn.close()
	return tt

def replace_last(source_string, replace_what, replace_with):
	head, sep, tail = source_string.rpartition(replace_what)
	return head + replace_with + tail

def to_float(in_str):
	return '{0:.4f}'.format(float(in_str)/1024 /1024 /1024)

def db_to_content(db , config_file):
	subject , content  = '' , ''
	content = 'Body = <html><head></head><body>'
	table = '<table border="1"> <tr>'
	for i in ['ID' , 'NAME' , '最终申请P' ,'最终申请VF', '最终申请HVMEM','使用IO' , '使用线程' , '使用最大内存', 'job_id' ,  'OLDP','OLDVF','OLDHVMEM' , 'MAXCYCLE', 'QSUBTIME','STATE']:
		table += '<th>{0}</th>  '.format(i)
	wrong_job = '<p>'
	flag = 0
	failed = 0 
	for i in query_all(db , 'SCRIPT'):
		table += '<tr>  '
		(ID,NAME,SCRIPT,QUEUE,P,VF,HVMEM,OTHER,OLDP,OLDVF,OLDHVMEM,MAXCYCLE,JOBID,QSUBTIME,FINISH,STATE) = i
		
		if P > OLDP or VF > OLDVF  :
			flag = 1 
			wrong_job += '&nbsp'.join(['资源申请过小','最终申请线程', P ,'原始申请线程', OLDP ,'最终申请内存', to_float(VF) , '原始申请内存' , to_float(OLDVF) , NAME]) + '<br>'
		if STATE != 'finished' : 
			flag = 1 
			failed = 1
			wrong_job += '&nbsp'.join(['任务反复投递失败' , NAME]) + '<br>'
		mm = query(db , 'JOB' , 'JOBID' , JOBID)
		_UIO , _UTHREAD , _UMAXVMEM , _JOBID = 0  , 0  , 0 , 0
		if len(mm) > 0 :
			[_JOBID , _NAME , _SCRIPT , _STATE , _P , _VF , _HVMEM , _OTHER , _NODE , _UCPU , _UMEM , _UIO , _UMAXVMEM, _UTHREAD, _START, _END , _FINISH , _KILL , _DIED] = mm[-1]
		
		if float(_UIO) > 10 or 0 < float(_UTHREAD) < float(OLDP)/2 or 0 < float(_UMAXVMEM) < float(OLDVF)/2 :
			flag = 1  
			wrong_job += '&nbsp'.join([ '资源申请过大' , 'IO:',_UIO ,'使用线程数', _UTHREAD ,'原始申请线程数',  OLDP , '使用最大内存', to_float(_UMAXVMEM) ,'原始申请内存', to_float(OLDVF) , NAME]) + '<br>'
		subject = replace_last( os.path.dirname(SCRIPT) , r'.qsub' , '')
		for i in [ID , NAME , P , to_float(VF), to_float(HVMEM),  str(_UIO) , str(_UTHREAD) , to_float(_UMAXVMEM) , str(_JOBID),  OLDP, to_float(OLDVF) , to_float(OLDHVMEM) , MAXCYCLE, QSUBTIME,STATE]: 
			table += '<td>{0}</td>  '.format(i)
		table += '</tr> '
	table += '</table>'
	wrong_job += '</p>'

	if flag == 1 : 
		content = content + table +  wrong_job +  '<p>{0}</p></body></html>'.format(subject)
		#print(wrong_job)
	else :
		content = content + table + '<p>{0}</p></body></html>'.format(subject)

	
	if failed : 
		subject = '[任务运行报告][失败]' + subject  
	else:
		subject = '[任务运行报告][成功]' + subject  
	email_file = '{0}.ini'.format(db)
	with open(email_file , 'w') as f_out :
		with open( config_file ) as f_in :
			for line in f_in :
				if line.startswith('Copy'):
					if flag == 1 : 
						f_out.write(line.rstrip() + ';taoliu@genome.cn\n')
					else :
						f_out.write('Copy =\n')
				else:
					f_out.write(line)
		f_out.write('[BODY]\n')
		f_out.write('Subject = {0}\n'.format(subject))
		f_out.write('Attachment = \n')
		f_out.write(content)
	os.chmod(email_file , 0o700)
	return email_file

def directory_full(db , config_file):
	subject , content , quota = '' , '' , ''
	for i in query_all(db , 'PARAMETER'):
		(ID  , MAXJOBS , INPUT , QUOTA )= i
		subject = INPUT 
		quota = QUOTA

	content = 'Body = <html><head></head><body><p>目录配额满了，请赶紧修改</p><p>{0} full ,quota is {1} ,  please to modify quota</p></body></html>'.format(subject , quota) 


	subject = '[错误：运行报告]' + subject
	email_file = '{0}.ini'.format(db)
	with open(email_file , 'w') as f_out :
		with open( config_file ) as f_in :
			for line in f_in :
				f_out.write(line)
		f_out.write('[BODY]\n')
		f_out.write('Subject = {0}\n'.format(subject))
		f_out.write('Attachment = \n')
		f_out.write(content)
	os.chmod(email_file , 0o700)
	return email_file


def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-full','--full',help='directory full ',dest='full', action='store_true')
	#parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	parser.add_argument('-m','--max_count' , help='max retry times',dest='retry', default = 0 , type=int)
	args=parser.parse_args()
	home_dir =  os.environ['HOME']
	if not os.path.isfile('{0}/.email/.email.txt'.format(home_dir)):
		print('''please config your email address first :
[DEFAULT]
Max_count = 5 #发送失败时重复发送次数
Sleep_time = 60 #发送失败时睡眠时间

[HEADER]
Addressor = **@genome.cm #邮件发送者
Password = **** #邮箱密码
Receiver = chengfangtu@annoroad.com; #接收者
Copy = yanxunsu@annoroad.com #抄送者
Server = smtp.exmail.qq.com #登陆服务器
Receive_server = pop.exmail.qq.com
''')
		sys.exit()
	else:
		if args.full:
			body_file = directory_full(args.input , '{0}/.email/.email.txt'.format(home_dir))
		else : 
			body_file = db_to_content(args.input , '{0}/.email/.email.txt'.format(home_dir))
		email = myemail.Email([body_file] , args.retry)
		email.send_email()

if __name__ == '__main__':
	main()
