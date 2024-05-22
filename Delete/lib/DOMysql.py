'''
连接lims数据库，从流程的config文件中获取到Lims数据库的连接信息，进行连接
执行单个功能的时候，按照table_name进行即可
'''
import os
import sys
import argparse
import mysql.connector
import re
pat1=re.compile('^\s+$')

__author__ = 'suyanxun'
__mail__ = 'yanxunsu@genome.cn'
pat1=re.compile('^\s+$')

class SQL():
    def __init__( self, config_file, port = 3306, db_type='local',charset='utf8' ):
        '''
        功能：连接和获取数据库连接
        参数：
        config_file: 配置文件，包含数据库连接信息
        port: 数据库端口，默认3306
        db_type: 数据库类型，默认local，可选['lims','local']
        charset: 数据库字符集，默认utf8
        '''
        self.config_file = config_file
        self.config_dic = read_config( self.config_file )
        if db_type == 'local':
            usr  = self.config_dic['sql_usr']
            pwd = self.config_dic['sql_pwd']
            port = self.config_dic['sql_port']
            host = self.config_dic['sql_host']
            database = self.config_dic['sql_db']
        else:
            usr  = self.config_dic['lims_sql_usr']
            pwd = self.config_dic['lims_sql_pwd']
            port = self.config_dic['lims_sql_port']
            host = self.config_dic['lims_sql_host']
            database = self.config_dic['lims_sql_db']
        try:
            self.cnx = mysql.connector.connect(user=usr, password=pwd, host=host, database=database, port = port, charset = charset)
            print( 'connect db {0}-{1} successed'.format(host, usr) )
            self.cursor = self.cnx.cursor()
        except mysql.connector.Error as err:
            print( 'connect db {0}-{1} failed'.format(host, usr) )
            print ( 'Error: {0}'.format( err ) )
            sys.exit(1)
        self.charset = charset
    
    def colname( self, table ):
        '''
        功能: 获取表列名
        参数：
            table: 表名
        返回值：
            name_list:list, 表列名
        '''
        cmd = 'describe {0}'.format( table )
        self.execute( cmd )
        name_list = [i[0] for i in self.cursor.fetchall() ]
        return name_list
    
    def select( self, table_name, col_list = '*', conditions = None ):
        '''
        功能: 查询表，默认返回所有值
        参数：
            col_list: list, 查询的列名
            conditions: list, 查询的条件,形如[ (colname1,value1), (colname2,value2), ... ]
        返回值:
            result: list, 查询结果,形如[("1","PM-"),("2","PM-2")]
        '''
        if col_list == '*':
            col = col_list
        else:
            col = ','.join( col_list )
        cmd = 'SELECT {0} from {1}'.format( col, table_name )
        if conditions:
            cmd = cmd + ' where {0} = "{1}"'.format( conditions[0][0], conditions[0][1] )
            if len( conditions ) > 1:
                for n,v in conditions[1:]:
                    cmd = cmd + ' and {0} = "{1}"'.format( n,v )
        self.execute( cmd )
        return self.cursor.fetchall()
    
    def insert( self, table_name, name_list, value_list ):
        '''
        功能: 插入数据
        参数：
            table_name: 表名
            col_list: list,要插入的列名list
            value_list: list,要插入的值,形如[(value1.1,value1.2,...), (value2.1,value2.2,...)]
        '''
        format_value = lambda x : '"{0}"'.format(str(x))
        cmd = 'INSERT INTO {0} ( {1} ) VALUES ({2});'.format(table_name , ",".join(name_list) , ",".join( map( format_value , value_list) ))
        self.execute(cmd)
    
    def update( self, table_name, value_list, conditions = None ):
        '''
        功能: 更新数据
        参数：
            table_name: 表名
            value_list: 要更新的值,形如[(colname1,value1), (colname2,value2), ... ]
            conditions: 更新条件: 形如[ (colname1,value1), (colname2,value2), ... ]
        使用方法：
        ziduan= []
        value = []
        value_list=list(zip(ziduan,value))
        update(table_name,, ziduan,value)
        '''
        cmd = 'UPDATE {0} SET '.format( table_name )
        for n,v in value_list:
            if n=='reads':
                cmd = cmd + '`{0}`="{1}", '.format( n,v )
            else :
                cmd = cmd + '{0}="{1}", '.format( n,v )
        cmd = cmd.rstrip( ', ' )
        if conditions:
            cmd = cmd + ' where {0} = "{1}"'.format( conditions[0][0], conditions[0][1] )
            if len( conditions ) > 1:
                for n,v in conditions[1:]:
                    cmd = cmd + ' and {0} = "{1}"'.format( n,v )
        self.execute( cmd )
    
    def delete( self, table_name, conditions ):
        '''
        conditions: 删除条件: 形如[ (colname1,value1), (colname2,value2), ... ]
        '''
        cmd = 'delete from {0} where {1}="{2}"'.format( table_name, conditions[0][0], conditions[0][1] )
        if len( conditions ) > 1:
            for n,v in conditions[1:]:
                cmd = cmd + ' and {0} = "{1}"'.format( n,v )
        self.execute( cmd )
    
    def execute( self, cmd, times = 3 ):
        if times > 0:
            try:
                self.cursor.execute(cmd)
                if cmd.startswith( ('INSERT' ,'UPDATE' , 'DELETE')):
                    self.cnx.commit()
            except mysql.connector.Error as err:
                print ('{0} 尝试倒数第{1}次失败'.format( cmd, times ) )
                print ( 'Error: {0}'.format( err ) )
                self.execute( cmd, times-1 )
        else:
            self.close()
            sys.exit()
    
    def close( self ):
        self.cnx.commit()
        self.cursor.close()
        self.cnx.close()
    
def read_config(config_file):
    '''
    功能：读取等号连接的config文件
    注意事项：
        1. 配置文件必须为等号连接的格式
        2. 配置文件中的key和value不能有空格
        3. 配置文件中的key和value不能有注释
        4. 配置文件中的key不能重复
    参数：
        config_file：配置文件，已key=value的格式配置
    返回值：
        dict:key=配置文件中等号前的信息，value=配置文件中等号后的信息
    '''
    config_dic = {}
    with open(config_file, 'r') as infile:
        for line in infile:
            if re.search(pat1, line ) or line.startswith('#') : continue
            tmp = line.rstrip().split('=',1)
            target = tmp[0].rstrip(' ')
            value = tmp[1].lstrip(' ')
            if target not in config_dic :
                config_dic[ target ] = value
            else :
                print("{0} is repeat in {1}".format(target, config_file))
    return config_dic

