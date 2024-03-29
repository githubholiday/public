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

class LIMS():
    def __init__( self, config_file, port = 3306, charset='utf8' ):
        self.config_file = config_file
        self.config_dic = self.read_config( )
        usr  = self.config_dic['sql_usr']
        pwd = self.config_dic['sql_pwd']
        port = self.config_dic['sql_port']
        host = self.config_dic['sql_host']
        database = self.config_dic['sql_db']
        try:
            self.cnx = mysql.connector.connect(user=usr, password=pwd, host=host, database=database, port = port, charset = charset)
            print( 'connect db {0}-{1} successed'.format(host, usr) )
            self.cursor = self.cnx.cursor()
        except mysql.connector.Error as err:
            print( 'connect db {0}-{1} failed'.format(host, usr) )
            print ( 'Error: {0}'.format( err ) )
            sys.exit(1)
        self.charset = charset
        
    def read_config(self):
        config_dic = {}
        with open(self.config_file, 'r') as infile:
            for line in infile:
                if re.search(pat1, line ) or line.startswith('#') : continue
                tmp = line.rstrip().split('=',1)
                target = tmp[0].rstrip(' ')
                value = tmp[1].lstrip(' ')
                if target not in config_dic :
                    config_dic[ target ] = value
                else :
                    print("{0} is repeat in {1}".format(target, self.config_file))
        return config_dic
    
    def colname( self, table ):
        cmd = 'describe {0}'.format( table )
        self.execute( cmd )
        name_list = [i[0] for i in self.cursor.fetchall() ]
        return name_list
    
    def select( self, table_name, col_list = '*', conditions = None ):
        '''
        默认返回所有值
        col_list: 查询的列名,list
        conditions: 条件,形如[ (colname1,value1), (colname2,value2), ... ]
        '''
        if col_list == '*':
            col = col_list
        else:
            col = ','.join( col_list )
        cmd = 'SELECT {0} from {1}'.format( col, table_name )
        if conditions:
            if "NULL" in conditions[0][2] :
                cmd = cmd + ' where {0} {1} {2}'.format( conditions[0][0], conditions[0][1], conditions[0][2] )
            else:
                cmd = cmd + ' where {0} {1} "{2}"'.format( conditions[0][0], conditions[0][1], conditions[0][2] )
            if len( conditions ) > 1:
                for n,r,v in conditions[1:]:
                    if "NULL" in v:
                        cmd = cmd + ' and {0} {1} {2}'.format( n,r,v )
                    else:
                        cmd = cmd + ' and {0} {1} "{2}"'.format( n,r,v )
        self.execute( cmd )
        return self.cursor.fetchall()
    
    #def insert( self, table_name, col_list, value_list ):
    def insert( self, table_name, name_list, value_list ):
        '''
        col_list: 要插入的列名list
        value_list: 要插入的值,形如[(value1.1,value1.2,...), (value2.1,value2.2,...)]
        '''
        format_value = lambda x : '"{0}"'.format(str(x))
        cmd = 'INSERT INTO {0} ( {1} ) VALUES ({2});'.format(table_name , ",".join(name_list) , ",".join( map( format_value , value_list) ))
        self.execute(cmd)
    
    def update( self, table_name, value_list, conditions = None ):
        '''
        value_list: 要更新的值,形如[(colname1,value1), (colname2,value2), ... ]
        conditions: 更新条件: 形如[ (colname1,value1), (colname2,value2), ... ]
        使用方法：
        ziduan=[]
        value = []
        value_list=list(zip(ziduan,value))
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
