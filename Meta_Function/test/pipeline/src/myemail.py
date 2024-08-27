'''
功能：根据config文件发送邮件
config格式：
[DEFAULT]
Max_count = 5 #发送失败时重复发送次数
Sleep_time = 60 #发送失败时睡眠时间

[HEADER]
Addressor = SCI_QC@annoroad.com #邮件发送者
Password = just4Use #邮箱密码
Receiver = chengfangtu@annoroad.com; #接收者
Copy = yanxunsu@annoroad.com #抄送者
Server = smtp.exmail.qq.com #登陆服务器
Receive_server = pop.exmail.qq.com

[PATHWAY]
Chugao_dir = /annoroad/data1/fqdata/tables/2016/sci/prefile/chugao
Shangji_dir = /annoroad/data1/fqdata/tables/2016/sci/prefile/shangji
FC_dir = /annoroad/data1/fqdata/tables/2016/sci/prefile/FC
BMS_dir = /annoroad/data1/fqdata/tables/2016/sci/prefile/BMS

[BODY]
Subject = 空邮件 #邮件主题
Attachment =  /annoroad/data1/bioinfo/PMO/suyanxun/test/test.txt #邮件附件
Body = 自动发邮件测试 #邮件内容，多行时请以空格开始每一行
 <加粗>
 不加粗

'''
#! /usr/bin/env python3
#coding: utf-8
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import utils, encoders
import mimetypes, sys,smtplib,socket,getopt
import argparse
import sys
import os
import re
import glob
import time
import configparser
import unittest

__author__='suyanxun'
__mail__= 'yanxunsu@annoroad.com'

pat1=re.compile('^\s+$')
#server = 'smtp.exmail.qq.com'
BLANK = '&nbsp'

class Email:
	def __init__(self,file_list , retry=0 ):
		self.file_list = file_list
		config = configparser.ConfigParser()
		for file in self.file_list:
			config.read(file , encoding='utf-8')
		self.addressor = config['HEADER']['Addressor']
		self.password = config['HEADER']['password']
		self.receiver = config['HEADER']['Receiver']
		self.copy = config['HEADER']['Copy']
		self.server = config['HEADER']['Server']
		self.subject = config['BODY']['Subject']
		self.body = self.get_body( config['BODY']['Body'] )
		if ( config['BODY']['Attachment'] != '' ):
			self.attachment = self.get_attachment( config['BODY']['Attachment'].encode('utf-8').decode() )
		else:
			self.attachment = ''
		if retry == 0:
			self.max_count = config['DEFAULT']['Max_count']
		else:
			self.max_count = retry
		self.sleep_time = config['DEFAULT']['Sleep_time']
	def get_body( self , body ):
		self.body = body
		#self.body = re.sub ( '<(.*?)>' , '<b>\\1</b>' , self.body )
		self.body = re.sub ( '\n' , '</br>' , self.body )
		return ( self.body )
	def get_attachment( self , file ):
		self.file = file.rstrip() 
		if ( not os.path.exists( self.file ) and self.file != '' ):
			print ('{0} not exists'.format( self.file ) )
			sys.exit(1)
		fd = open( '{0}'.format( self.file ) , 'rb' )
		mimetype,mimeencoding=mimetypes.guess_type( self.file.split( '/' )[-1] )
		if ( mimeencoding is None ) or ( mimetype is None ):
			mimetype = 'application/octet-stream'
			maintype,subtype = mimetype.split('/')
		if maintype=='text':
			attachment=MIMEText( fd.read(), _subtype=subtype, _charset='utf-8' )
		else:
			attachment = MIMEBase( maintype , subtype )
			attachment.set_payload( fd.read() )
			encoders.encode_base64( attachment )
			attachment.add_header( 'Content-Disposition' , 'attachment' , filename = self.file.split('/')[-1] )
		fd.close()
		return attachment
	def msginfo( self ):
		msg=MIMEMultipart()
		msg['From'] = self.addressor
		msg['To'] = self.receiver
		if ( self.copy != '' ):
			msg['Cc'] = self.copy
		msg['Date'] = utils.formatdate(localtime=1)
		msg['Message-ID'] = utils.make_msgid()
		msg['Subject'] = self.subject
		print(self.body)
		self.body = MIMEText( self.body , 'html' , 'utf-8' )
		msg.attach( self.body )
		if (self.attachment != ''):
			msg.attach( self.attachment )
		return msg.as_string()
	def send_email( self ):
		count = 0
		flog = True
		smtp=smtplib.SMTP(self.server)
		while ( count < int( self.max_count ) ):
			try:
				flog = True
				try:
					flog = True
					smtp.login( self.addressor , self.password )
				except Exception as e:
					print( e )
					print( 'Email login failed:' )
					flog=False
				if ( self.copy != ''):
					receiver=self.receiver+';'+self.copy
				else :
					receiver=self.receiver
				smtp.sendmail( self.addressor , receiver.split(';') , self.msginfo())
			except Exception as e:
				print (e)
				print ('Send mail failed')
				flog=False
			else:
				print ("Send to {0} receiver sucessendly".format(len(receiver)))
			if ( flog ):
				break
			else:
				count += 1
				if ( count == int( self.max_count ) ):
					print("Try {0} and fail").format( self.max_count )
					sys.exit()
				time.sleep( int( self.sleep_time ) )
class TestEmail(unittest.TestCase):
	def test_email(self):
		config = ['/annoroad/data1/bioinfo/PMO/suyanxun/test/example1.ini']
		a = Email( config )
		a_subject = a.subject
		#a_attachment = a.get_attachment( a_attachment_file )
		self.assertTrue( a_subject , '自动发邮件测试' )
		#self.assertTrue( a_body , 'Hello: </br>{0}<b>dh</b><br/>{0}fv<br/>'.format ( BLANK*6 ) )
		#self.assertTrue( a_attachment_file , config )
		a.send_email()

if __name__ == '__main__':
	   unittest.main()
