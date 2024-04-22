##### cat example.py 
#! /usr/bin/env python3
import argparse
import sys
import os
import re
import datetime
import glob
import time
bindir = os.path.abspath(os.path.dirname(__file__))
filename=os.path.basename(__file__)

__author__='tu chengfang '
__mail__= 'chengfangtu@genome.cn'


class Log():
	def __init__( self, filename, funcname = '' ):
		self.filename = filename 
		self.funcname = funcname
	def format( self, level, message ) :
		date_now = datetime.datetime.now().strftime('%Y%m%d %H:%M:%S')
		formatter = ''
		if self.funcname == '' :
			formatter = '\n{0} - {1} - {2} - {3} \n'.format( date_now, self.filename, level, message )
		else :
			
			formatter = '\n{0} - {1} - {2} -  {3} - {4}\n'.format( date_now, self.filename, self.funcname, level, message )
		return formatter
	def info( self, message ):
		formatter = self.format( 'INFO', message )
		sys.stdout.write( formatter )
	def debug( self, message ) :
		formatter = self.format( 'DEBUG', message )
		sys.stdout.write( formatter )
	def warning( self, message ) :
		formatter = self.format( 'WARNING', message )
		sys.stdout.write( formatter )
	def error( self, message ) :
		formatter = self.format( 'ERROR', message )
		sys.stderr.write( formatter )
	def critical( self, message ) :
		formatter = self.format( 'CRITICAL', message )
		sys.stderr.write( formatter )
def get_value(in_str, target ):
	if target in in_str:
		tmp = in_str.split('"')
		value = tmp[1]
		return value
	else:
		return 0

def read_email_feedback_txt(email_feedback_txt):
	with open(email_feedback_txt) as f:
		for line in f:
			lines = line.strip()
			if line.startswith('finish'):
				line = line.strip()
				tmp = line.split('=')[1]
				if tmp == 'YES':
					return True
				else:
					return False

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--indir',help='input file',dest='indir',required=True)
	parser.add_argument('-o','--outfile',help='input file',dest='outfile',required=True)
	args=parser.parse_args()
	
	output = open(args.outfile,'w')
	filter_result_dir_pattern = '{0}/*/*/Filter_Result/'.format(args.indir )
	filter_result_dir_list = glob.glob( filter_result_dir_pattern )
	for filter_result_dir in filter_result_dir_list:
		email_feedback_txt = '{0}/email/email_feedback.txt'.format( filter_result_dir )
		data_finish = read_email_feedback_txt(email_feedback_txt)
		print(email_feedback_txt)
		if data_finish == True:
			output.write( filter_result_dir + '\n' )
			file_create_time = time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getctime(email_feedback_txt)))
			output.write( '\t'.join([filter_result_dir, file_create_time]) + '\n' )
	output.close()

	
if __name__ == '__main__':
	main()