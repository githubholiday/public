#! /usr/bin/env python3
import configparser
import datetime
import os
import sys
import unittest
import re
filename = os.path.basename(__file__)
pat1=re.compile('^\s+$')
__author__='Tu Chengfang'
__mail__= 'chengfangtu@genome.cn'
pat=re.compile('^\s$')
pat2 = re.compile('\[(\S+)\]')
date_now = datetime.datetime.now().strftime('%Y%m%d %H:%M:%S')

class Config():
	def __init__(self, config_file):
		self.config_file = config_file
		self.check_file_exist()

	def check_file_exist(self):
		if not os.path.exists(self.config_file ):
			error_message = '{0}文件中不存在'.format(self.config_file)
			sys.stderr.write('\n{0} - {1} - {2} - {3} \n'.format( date_now, filename, 'ERROR', error_message ))
			sys.exit(1)
			
	def read_config( self ):
		config_dic = {}
		with open( self.config_file ) as infile :
			for line in infile :
				if re.search(pat1, line ) or line.startswith('#') : continue
				tmp = line.rstrip().split('=')
				target = tmp[0].rstrip(' ')
				value = tmp[1].lstrip(' ')
				if target not in config_dic :
					config_dic[ target ] = value
				else :
					error_message = '{0}文件中的target有重复:{1}'.format(self.config_file, target)
					sys.stderr.write('\n{0} - {1} - {2} - {3} \n'.format( date_now, filename, 'ERROR', error_message ))
					sys.exit(1)
		return config_dic

	def get_value( self, target ):
		config_dic = self.read_config()
		try :
			value = config_dic[target]
			return value
		except :
			error_message = '{0} not in config,Please check'.format( target )
			sys.stderr.write('\n{0} - {1} - {2} - {3} \n'.format( date_now, filename, 'ERROR', error_message ))
			sys.exit(1)

class myconf(configparser.ConfigParser):
	def __init__(self,defaults=None):
		configparser.ConfigParser.__init__(self,defaults=None,allow_no_value=True)
	def optionxform(self, optionstr):
		return optionstr


class Config_Parser( ) :
	def __init__( self, config_file ) :
		self.config_file = config_file 
		self.config = self.read_config()
		
	def read_config( self ) :
		config = myconf()
		config.read( self.config_file )
		return config

	def all_block( self, title, head ) :
		return self.config[title][head]
		
	def return_block_list( self, title ) :
		try :
			data = self.config[ title]
		except :
			return []
		for rn in data :
			info = data[rn].rstrip().split('\t')
			yield info

def read_analysisConf( config ) :
    dict = {}
    header = ''
    with open ( config, 'r' ) as IN :
        for line in IN :
            if line.startswith('#') or re.search( pat, line ) : continue
            if line.startswith('[') :
                match = pat2.search(line)
                if match :
                    header = match.group(1)
                    dict[header] = {}
            else :
                if header == 'Para':
                    key, value = line.replace(' ', '').rstrip().split('=',1)
                    dict[header][key] = value
                else :
                    list = line.rstrip().split('\t')
                    num = len( dict[header] )
                    dict[header][num] = list
    return dict