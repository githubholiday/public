import argparse
import sys
import re
import os
from datetime import datetime
from dateutil import rrule
import time
import glob
import configparser
from collections import defaultdict, OrderedDict
bindir = os.path.abspath(os.path.dirname(__file__))
filename = os.path.basename(__file__)
sys.path.append('{0}/../../lib'.format(bindir))
from readconfig import Config
import logger
my_log = logger.Log(filename)
# JUMP_SCRIPT = '{0}/../jump_email/jump_send_email.py'.format( bindir )
__author__='Yang Zhang'
__mail__= 'yangzhang@genome.cn'
__doc__='''
功能：
这是一个生成邮件内容的函数。
使用方法：
from prepareEmail import Email
email_prepare = Email(email_ini, email_file , person, mail)
email_dir, email_ini_out = email_prepare.email_info_write()
然后就可以发送邮件了：
send_qc_email_cmd = '{PYTHON} {bindir}/send_email.py -i {email_ini_out}'.format(PYTHON=PYTHON, bindir=bindir, email_ini_out=email_ini_out )

参数：
email_ini:发送邮件的参数，详见config/pipe_email.ini
email_file:输入的表格，用于后续编译为发送邮件的内容
mail:收件人，项管名字下的所有邮箱"a@genome.cn;b@genome.cn"
'''

def std( level, message ):
    now_time = time.strftime("%Y-%m-%d %H:%M:%S")
    string = '{0} - {1} - {2} - {3}\n'.format( now_time, filename, level, message )
    if level == 'ERROR':
        sys.stderr.write( string )
    else:
        sys.stdout.write( string )

def myrun( cmd, prompt ='' ) :
    if not prompt :
        prompt = cmd
    if os.system( cmd ) == 0 :
        if 'mkdir' in cmd :
            time.sleep(1)
        my_log.info( '{0} ## run sucessfully!'.format( prompt ))
        return 1
    else :
        my_log.error( '{0} ## run failed!'.format( prompt ))
        return 0
    
class Email():
    def __init__( self, email_ini , email_file, person):
        self.email_ini = email_ini
        self.email_file =  '{0}'.format( email_file.replace("xls","ini"))
        self.email_dir = os.path.dirname(email_file)
        self.person = person
        self.day = datetime.now().strftime('%Y-%m-%d')
        self.target = 'Filter_Delete'
        self.file = email_file
        # self.get_dir()
        self.mkdir()
        
    def config_read( self ) :
        self.pipe_email_conf = configparser.ConfigParser()
        self.pipe_email_conf.read( self.email_ini )
        return self.pipe_email_conf
          
    def mkdir( self ) :
        self.scan_email_dir = '{self.email_dir}'.format( self=self )
        if not os.path.exists( self.scan_email_dir ) :
            mkdir_cmd = 'mkdir -p {0}'.format( self.scan_email_dir )
            run_state = myrun( mkdir_cmd )
            if not run_state :
                my_log.error('mkdir {0} fail'.format(self.scan_email_dir))
                sys.exit(1)
        
    def get_email_setting( self ) :
        email_config = self.config_read()
        section = self.target
        self.max_count = email_config['DEFAULT']['Max_count']
        self.sleep_time = email_config['DEFAULT']['Sleep_time']
        self.addressor = email_config['HEADER']['Addressor']
        self.password = email_config['HEADER']['Password']
        self.server = email_config['HEADER']['Server']
        self.receive_server = email_config['HEADER']['Receive_server']
        self.subject = email_config[section]['Subject']
        self.receiver = email_config[section]['Receiver']
        # self.receiver = self.mail
        self.copy = email_config[section]['Copy']

    def get_email_body( self ):
        body_file =  '{0}/email.txt'.format( self.email_dir ) 
        # print(body_file)
        if not body_file :
            my_log.error('The file with the warning deletion message is not in path:{0}'.format(self.email_dir ))
        else:
            self.subject = '{self.subject}'.format( self=self )
        self.email_body = self.body( body_file )
    
    def body( self, body_file ):
        email_body = '<html><head></head><body><P>{0} ：</P><P>~~</P></P><P>{1}</P>'.format(self.person, "请关注附件中的项目，核对情况，及时交付，本地路径结果即将删除")
        # email_file_dic = OrderedDict()
        email_file_dic = {'body_file':[body_file,'项目统计']}  
        file_list = ['body_file']
        # print(email_file_dic)
        for target in file_list :
            file_path = email_file_dic[target][0]
            # print(file_path)
            # email_title = email_file_dic[target][1]
            if os.path.exists( file_path ):
                # email_body += '<p><span>{0}:</span></p>'.format(email_title)
                with open( file_path, 'r' ) as infile:
                    email_body += '<table style="BORDER-COLLAPSE: collapse" borderColor=#000000 height=40 width=1250 cellPadding=1 align=center border=1px">'
                    for line in infile:
                        #print(line)
                        email_body += '<tr>'
                        temp = []
                        tmp = line.rstrip().split('\t')
                        for j in tmp :
                            j = j.replace('(%)','') # 发邮件的时候如果有%会报错
                        email_body += ''.join( temp )
                        email_body += '</tr>'
                    email_body += '</table>'
            else:
                std( 'WARNNING', ' {0} 不存在'.format( file_path ) )
        return email_body

    def email_info_write( self,now_time ):
        self.get_email_setting()
        self.get_email_body()
        self.email_ini_out = self.email_file
        with open(self.email_ini_out, 'w') as outfile :
            outfile.write('''
[DEFAULT]
Max_count={self.max_count}
Sleep_time={self.sleep_time}
[HEADER]
Addressor = {self.addressor}
Password = {self.password}
Receiver  = {self.receiver}
Copy = {self.copy}
Server = {self.server}
Receive_server = {self.receive_server}
[BODY]
Subject={self.subject}-{now_date}
Attachment={self.file}
Body={self.email_body}                        
            '''.format( self=self,now_date=now_time ))
        return self.email_dir, self.email_ini_out
