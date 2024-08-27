#! /usr/bin/env python3
# -*- coding: utf-8 -*-  
import argparse
import sys
import os
import re
import paramiko
import sqlite3
import pexpect
bindir = os.path.abspath(os.path.dirname(__file__))
import myemail

__author__='Liu Tao'
__mail__= 'taoliu@annoroad.com'

pat1=re.compile('^\s+$')
IP='192.168.60.188'
IP2 = '192.168.13.20'
PORT='50732'
PYTHON='/usr/local/bin/python3'
SCRIPT='/home/sci-qc/bin/qsub_sge/send_email.py'

def check_jump_success(user):
	cmd = "ssh {0}@{1} -p {2}".format( user , IP , PORT )
	print(cmd)
	child = pexpect.spawn( cmd )
	state = [True , True , True , False , False]
	index = child.expect([ "Last login:", "Welcome", pexpect.EOF, "password", pexpect.TIMEOUT ], timeout = 30)
	if index < len(state):
		result = state[index]
	else:
		result = False
	child.close( force=True )
	return(result)

def exists_email_config(user):
	cmd = "ssh {0}@{1} -p {2}".format( user , IP , PORT )
	print(cmd)
	child = pexpect.spawn( cmd )
	state = [True , True , True , False , False]
	index = child.expect(["Welcome", "Last login:", "Welcome", pexpect.EOF, "password", pexpect.TIMEOUT ], timeout = 30)
	print(index)
	if index < len(state):
		result = state[index]
	else:
		result = False
	cmd = "test -f ~/.email/.email.txt && echo ok"
	print(cmd)
	tt = mysendline( child , cmd )
	child.close( force=True )
	return(tt)

def mysendline(child , cmd ) :
	print(cmd)
	child.sendline(cmd)
	index =  child.expect( ['tt' , 'ok' , pexpect.EOF , pexpect.TIMEOUT ] , timeout = 150)
	print(index)
	if index == 0 :
		return True
	else :
		sys.exit("{0} failed".format(cmd))


def scp_db_file( user , old_db , full_flag):
	cmd = "ssh {0}@{1} -p {2}".format( user , IP , PORT )
	child = pexpect.spawn( cmd )
	mysendline( child , "mkdir -p ~/$$ && echo ok" )
	mysendline( child ,  "scp {0}@{1}:{2} ~/$$/  && echo ok ".format(user , IP2, old_db))
	if full_flag:
		mysendline( child  , "{0}  {1} -i ~/$$/{2} -f   && echo ok ".format(PYTHON , SCRIPT , os.path.basename(old_db)))
	else:
		mysendline( child ,  "{0}  {1} -i ~/$$/{2}  && echo ok ".format(PYTHON , SCRIPT , os.path.basename(old_db)))
	child.close( force=True )

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-full','--full',help='directory full ',dest='full', action='store_true')
	#parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	#parser.add_argument('-m','--mm',help='output file',dest='mm',action='store_false')
	args=parser.parse_args()
	#user = getpass.getuser()
	user = 'sci-qc'
	if not check_jump_success(user):
		sys.exit("免密跳转失败，请配置")
	if not exists_email_config(user):
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
	scp_db_file(user , os.path.abspath(args.input) , args.full)

if __name__ == '__main__':
	main()
