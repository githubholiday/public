import os
import sys
import datetime
import configparser

def color_print( content, backgroud='40', font_color= '31', asc_control = '0' ) :
    '''
    use print with color
    backgroud : 背景色[40-49]
    font_color: 字体颜色[30-39]
    asc_control: 控制码[0-8]
    默认的为黑色背景，字体颜色为红色，并关闭所有属性
    backgroud= range(40,50)
    font_color = range(30,40)
    for b in backgroud:
        for f in font_color:
            print("b:{0} f:{1}".format(b,f))
            print("\033[{0};{1}m{2}\033[{3}m".format( b, f, "测试", '0' ))
    '''
    print("\033[{0};{1}m{2}\033[{3}m".format( backgroud, font_color, content, asc_control ))
    
def myinput( content, empty="no" ) :
    '''
    通过手动输入信息
    content:提示信息
    empty:是否可以为空，默认是不可以，必须填入信息
    '''
    result = ''
    while (1) :
        try :
            result = input ('\033[40;35m>>>{0}\033[0m\n'.format( content  ))
            result = result.rstrip()
            if not result and empty=="no":
                color_print( '请重新输入')
            else :
                if 'exit' in result :
                    sys.exit(0)
                    color_print( '程序退出')
                else :
                    return result
        except SystemExit :
            color_print('程序退出')
            sys.exit(0)
        except :
            color_print('请重新输入')
            
def time_stamp():
    '''
    返回时间戳，以ms为单位，可以使用时间戳转化工具查看具体时间
    '''
    import time 
    now_time = int(time.time()*1000)
    return now_time
    
def uuid():
    '''
    用于生成uuid
    '''
    import uuid
    uuid = str(uuid.uuid1()).replace('-','')
    return uuid
    
    
# class Read_Config(configParser) :
    # def __init__( self, config_file ):
        # configParser.ConfigParser.__init__(self,defaults=None,allow_no_value=True)
        # self.config_file = config_file
        # self.config = self.read_config( )
class Read_Config() :
    def __init__( self, config_file ):
        self.config_file = config_file
        self.config = self.read_config( )
    
    def read_config( self ) :
        config = configparser.ConfigParser()
        config.read( self.config_file)
        return config
        
class ReadReleaseConf():
    def __init__(self, release_conf):
        self.release_conf = release_conf
        release_conf_object = Read_Config(self.release_conf )
        self.config = release_conf_object.config
    def print_type( self ):
        '''
        retrun ：
        type_print_dic = {'1':'UBUNTU'}
        '''
        type_print_dic= {}
        for num in self.config['TYPE'] :
            type = self.config['TYPE'][num]
            if num not in type_print_dic:
                type_print_dic[num] = type
        return type_print_dic

    def select_info( self , key) :
        '''
        key:UBUNTU,CLOUD_DELIVERY 等
        ubuntu_user_dic:0=安诺云   annoyun 10.1.4.94   123 oss://annoroad-cloud-product/user/  10.4.1.94
        i=tmp[0]    tmp[1]  tmp[2]  tmp[3]  tmp[4]  tmp[5]
        info_print=['[1]xiazai1']
        info_user_dic={'xiazai1':{'name':'xiazai1','user':'root','ip':'192.168.','passwd':'123465','dir':'/media','node':'c008,c0019'}}
        num_user_dic={'1':'xiazai1'}
        '''
        key = key.upper()
        info_print = []
        num_user_dic = {}
        info_user_dic = {}
        info_num = self.config.options( key)#得到等号前的内容
        #print(info_num)
        for i in info_num :
            info_list = self.config[ key ][i] #得到等号后面的内容
            tmp = info_list.split('\t')
            name = tmp[0]
            out = '[{0}]{1}'.format( i, name)
            if out not in info_print :
                info_print.append( out )
            if name not in info_user_dic :
                info_user_dic[ name ] = {}
            info_user_dic[ name ]['name'] = name
            info_user_dic[ name ]['user'] = tmp[1]
            info_user_dic[ name ]['ip'] = tmp[2]
            info_user_dic[ name ]['passwd'] = tmp[3]
            info_user_dic[ name ]['dir'] = tmp[4]
            info_user_dic[ name ]['node'] = tmp[5]
            info_user_dic[ name ]['type'] = tmp[6]
            if i not in num_user_dic:
                num_user_dic[i] = name
        return info_print , num_user_dic, info_user_dic
    
    def select_other_info( self , key) : #except ubuntu
        other_dic = {}
        key = key.lower()
        tmp = self.config[ 'OTHER'][ key ].split('\t')
        if key not in other_dic :
            other_dic[ key ] = {}
        other_dic[ key ]['ip'] = tmp[0]
        other_dic[ key ]['dir'] = tmp[1]
        other_dic[ key ]['node'] = tmp[2]
        return other_dic
    def get_cmd( self, key, target ):
        key = key.lower()
        cmd = self.config[ key ][target] #如annoyun,upload
        return cmd

