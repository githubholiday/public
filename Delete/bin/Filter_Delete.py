#!/usr/bin/env python3
import os
import sys
import re
import argparse
import random
import glob
import datetime
# from datetime import datetime
from dateutil import rrule
import time
import configparser
bindir = os.path.abspath(os.path.dirname(__file__))
sys.path.append( '{0}/../lib'.format( bindir ))
from prepareEmail import Email
# from Lims_SQL import LIMS
import DOMysql
from readconfig import Config
from robot import Robots
import logger
filename = os.path.basename(__file__)
my_log = logger.Log(filename)
# from PipMethod import myconf,mkdir,generateShell
# from pathlib import Path


__author__='Yang Zhang'
__mail__= 'yangzhang@genome.cn'
__doc__='''
功能：
1、扫描数据库tb_deletion_info中数据交付是否完成，更新状态；（scan）
    提取数据库中删除类型为 “交付删除”，交付是否完成为 0（否）的项目，进行以下判断：
    （1）如果是hw（华为云），检查交付目录下是否有Finish.log文件，有的话，则更新交付完成时间为当前时间，交付是否完成为1（完成）
    （2）如果是ali（阿里云），检查交付目录下是否有Finish.log文件，有的话，则更新交付完成时间为当前时间，交付是否完成为1（完成）
    （3）如果是others（硬盘类），检查交付开始时间，如果已经>15天，则更新交付完成时间为当前时间，交付是否完成为1（完成）
    * 交付删除的含义是：已经发起过交付的命令。预期就是正在交付或者已经交付成功的项目。
    * 超期删除的含义是：没有发起过交付，但是已经超过了本地保留时间的项目。

2、根据配置文件中的时间（21天）进行超期提醒，把所有项目发送到项管的群邮箱 (email)，附件形式发送，需要转码。
    （1）扫描本地的Analysis_*目录，将扫到的目录存在outdir/all_local_dir_file.txt文件中
    （2）根据（1）中的文件，逐个目录判断：获得项目号，创建时间，删除时间（创建时间+30天，30可以在config中修改），超期判断；
        超期判断：
        1. 大前提: 本地目录创建时间距离当前时间超过21天
        2. 然后根据中间表信息进行判断
            1）云上无该项目记录：记录项目为超期删除，并邮件提醒
            2）云上有该项目，之前的记录中“是否删除”均为是，则记录项目为超期删除，并邮件提醒
            3）云上有该项目，之前的记录中“是否删除”含有否，只要有一个删除类型为 交付删除，直接跳过
            4）云上有该项目，之前的记录中“是否删除”含有否，无交付删除，邮件提醒
    （3）生成超期提醒的项目列表，根据项管名字分成多个文件写到outdir/email.person.txt
    （4）根据（3）中的文件，生成邮件正文outdir/email.person.ini
    （5）发送邮件，如果成功生成outdir/email.person.ini.sign标志，检查标志，如果不存在则重复发送3次，如果最后还是失败，就发送机器人提醒
    （6）【过滤项目数据超期删除提醒失败！！！失败的是：】项管 邮件发送命令


使用方法：（238集群）
/annoroad/data1/software/bin/miniconda/envs/python3_base/bin/python3 Filter_Delete.py -o outdir

参数说明：
-o outdir 为输出目录 ， 可以直接完成功能1和2。
-t True ， 进行测试 ， 测试与正式不同的在于这两个文件都是测试版本：
            pipe_email_test.ini（收件人信息）；
            tb_deletion_info_test（要修改信息的数据库）
-e False ， 只执行功能1 
-s False ， 只执行功能2 
其他相关参数暴露于../config/config.txt中。
'''

def myrun( cmd, prompt ='' ) :
    ''' 
    功能: 反馈linux命令执行结果
    返回值: 1是成功, 0是失败
    注意事项：无论成功与否, 都不会中断
    '''
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
    
def mkdir( dir ) :
    '''
    功能：创建目录
    返回值：无
    注意事项：创建目录失败会中断退出
    '''
    if not os.path.exists( dir ) :
        mkdir_cmd = 'mkdir -p {0}'.format( dir )
        run_state = myrun( mkdir_cmd )
        if not run_state :
            my_log.error('mkdir {0} fail'.format(dir))
            sys.exit(1)

def InterTime(Time,TimeFormat=''):
    ''' 
    功能：计算时间间隔。
    返回值：返回现在-输入时间之间的时间间隔，值为整数的天数
    注意事项：
    1、输入时间需要早于现在，否则返回0
    2、如果使用的Time是datetime格式，则不需要给TimeFormat（默认）；
       如果是其他格式，则需要通过TimeFormat参数给对应的格式形式
    '''
    if TimeFormat:
        try:
            start_time = datetime.datetime.strptime(Time,TimeFormat)
        except:
            start_time = datetime.datetime.strptime(Time,'%Y-%m-%d %H-%M-%S')
    else:
        start_time = Time
    now_time = datetime.datetime.now()
    days = rrule.rrule(freq=rrule.DAILY, dtstart=start_time, until=now_time).count()
    return days

def GetStatus(account,cloud_address,OBS_CMD,OSS_CMD):
    ''' 
    功能：查看华为云/阿里云上数据交付是否完成 
    返回值：True/False
    注意事项：可能网络原因或者config文件失效 导致扫描失败
    参数说明：
    account:hw,ali
    cloud_address:云上地址，该地址下有完成标志文件 Finish.log
    OBS_CMD: 执行的华为云命令
    OSS_CMD: 执行的阿里云命令
    '''
    if cloud_address[-1] == '/':
        Dir = cloud_address + "Finish.log"
    else:
        Dir = cloud_address + "/Finish.log"
    # print(Dir)
    number = 0
    if account == 'hw':
        p=os.popen("{0} {1}".format(OBS_CMD, Dir))
        # print("{0} {1}".format(OBS_CMD, Dir))
        for line in p.readlines():
            # print(line)
            if line.startswith("File number is:"):
                number = int(line.split(":")[1].strip())
            else:continue        
    elif account == 'ali':
        p=os.popen("{0} {1}".format(OSS_CMD, Dir))
        # print("{0} {1}".format(OSS_CMD, Dir))
        for line in p.readlines():
            # print(line)
            if line.startswith("Object Number is: "):
                number = int(line.split(":")[1].strip())
            else:continue
    if number == 1:
        FinishLog=True
    else:
        FinishLog=False
    return FinishLog # True/False

def Scan_Status(lims_config_file,table,config,config_dir,TimeFormat,ReleaseDay):
    ''' 
    功能：扫描中间表 tb_deletion_info 并标记交付成功与否的状态和时间 
    返回值：无，直接在中间表进行修改
    注意事项：如果没有更新数据库，会将状态记录在日志中，需要排查是否有问题
    参数说明：
    lims_config_file: 中间表的账号密码
    table: 中间表名字，这里用的是 tb_deletion_info
    config: 获取云上文件扫描命令的配置文件
    config_dir: 上述config所在路径
    TimeFormat: 时间格式化形式
    ReleaseDay: 15天，判断硬盘交付是否完成的时间长度
    '''
    lims_db_do = DOMysql.SQL( lims_config_file )
    lims_record = lims_db_do.select(table,conditions=[('delete_type', "交付删除"),('delivery_bool',0)]) # delivery_bool为0,1，整数，但是筛选的时候可以是0，也可以是'0'，么有影响
    OBS_CMD = config.get_value('obs_cmd').replace('bindir', config_dir)
    OSS_CMD = config.get_value('oss_cmd').replace('bindir', config_dir)

    for i in lims_record:
        account = i[6]
        cloud_address = i[3]
        start_time = i[2]
        project_id = i[1]
        FinishLog = GetStatus(account,cloud_address,OBS_CMD,OSS_CMD)
        days = InterTime(start_time,TimeFormat)
        now_time = datetime.datetime.now()
        UpdateFlag = False

        if account == 'others' and days > ReleaseDay: 
            UpdateFlag = True
        else:
            if account == 'ali' and FinishLog:
                UpdateFlag = True
            elif account == 'hw' and FinishLog:
                UpdateFlag = True
            else:
                my_log.info("未更新数据库：交付项目为 {0},\n云上路径为 {1}, \n交付账户为{2}, \n交付时间为{3}, \n云上Finish.log是否存在{4}\n".format(project_id, cloud_address, account, now_time.strftime(TimeFormat), FinishLog ))
        if UpdateFlag:
            delivery_bool = 1 # 是
            finish_time = now_time.strftime(TimeFormat)
            lims_db_do.update(table, [('delivery_bool',delivery_bool),('delivery_finish_time',finish_time)], [('project_id',project_id),('cloud_address',cloud_address)])

def Get_local_dir(CMD, outfile):
    '''
    功能：扫描本地过滤路径(Analysis_*）
    返回值：存储扫到的路径到输出文件
    注意事项：无
    '''
    p = os.popen(CMD) 
    OUT = open(outfile , 'w')
    for line in p.readlines():
        if line == '^$':continue
        else:
            OUT.write(line+"\n")
    OUT.close()

def Scan_local_dir(WarningDay,lims_config_file,table,local_dir_file,TimeFormat,DeleteDay):
    '''
    功能：根据local_dir_file文件，逐个判断是否有超期项目
    返回值：需要提醒的项目列表，根据删除时间排序，升序 [项目号，项目名称，项管名字，项管邮箱，本地目录创建时间，删除时间，本地目录（具体到Analysis_*）]
    注意事项：
    1. 项管的信息来自 tb_info_sequence_bill ，根据项目号进行搜索
    2. 如果出现异常判断，会记录为unknown
    超期判断标准：
    1. 大前提: 本地目录创建时间距离当前时间超过21天
    2. 然后根据中间表信息进行判断
        1）云上无该项目记录：记录项目为超期删除，并邮件提醒
        2）云上有该项目，之前的记录中“是否删除”均为是，则记录项目为超期删除，并邮件提醒
        3）云上有该项目，之前的记录中“是否删除”含有否，只要有一个删除类型为 交付删除，直接跳过
        4）云上有该项目，之前的记录中“是否删除”含有否，无交付删除，邮件提醒

    ''' 
    local_sql = DOMysql.SQL( lims_config_file )
    lims_db = DOMysql.SQL( lims_config_file, db_type='lims' )
    now_time = datetime.datetime.now().strftime(TimeFormat)
    name_list = ['project_id','delete_type','delivery_start_time','delivery_bool','delivery_finish_time','analysis_bool','local_delete_bool','create_time']
    p = open(local_dir_file,'r')
    delete_list = []
    for line in p:
        local_dir = line.strip()
        # print(local_dir)
        if not local_dir: # 防止空行 
            continue
        else:
            if not os.path.exists(local_dir):  continue # 防止处理到这步的时候，路径已经被删除了
            create_time = datetime.datetime.fromtimestamp(os.stat(local_dir).st_mtime) # datetime格式(2024,3,14,2,4,5)
            days = InterTime(create_time)
            delete_time = (create_time + datetime.timedelta(days=int(DeleteDay))).strftime(TimeFormat) 
            create_time = create_time.strftime(TimeFormat)
            lists = local_dir.split("/")
            project = [i for i in lists if "PM-" in i]
            if len(project) == 1 :
                projects = project[0].split("_")
                projects = [i for i in projects if "PM-" in i]
                project_ID = projects[0]
            elif len(project) == 0:
                print("这个目录下没有以PM-开头的层级，不能确认项目号，请核对："+local_dir)  # 以unknown记录异常情况
                project_ID = "unknown"
            else :
                print("这个目录下不止一个PM-开头的层级，不能确认项目号，请核对："+local_dir)
                project_ID = "unknown"
            try:
                record = lims_db.select("tb_info_sequence_bill", col_list=['latest_project_user_name','latest_project_user_email','task_name'], conditions=[('project_code',project_ID)])
                project_user_name = record[0][0]
                project_user_mail = record[0][1]
                project_name = record[0][2]
            except:
                project_user_name = "unknown"
                project_user_mail = "-"
                project_name = "-"
            value_list = [project_ID,'超期删除',create_time,1,now_time,0,0,now_time]
            warn_list = [project_ID,project_name,project_user_name,project_user_mail,create_time,delete_time,local_dir]
                # [项目号，项目名称，项管名字，项管邮箱，本地目录创建时间，删除时间，本地目录（具体到Analysis_*）]
            if days <= WarningDay :
                continue
            else :
                lims_record = local_sql.select(table,conditions=[('project_id', project_ID)])
                if len(lims_record) > 0: # 判断项目是否存在
                    local_delete_bool_list = [i[10] for i in lims_record]
                    if 0 in local_delete_bool_list:  # 是否删除为0,1，类型为整数
                        delete_type_list = [i[5] for i in lims_record if i[10] == 0]  
                        if '交付删除' in delete_type_list:
                            continue
                        else:
                            delete_list.append(warn_list)
                    else:
                        local_sql.insert(table,name_list,value_list)
                        delete_list.append(warn_list)
    p.close()
    # delete_list_sort = sorted(delete_list,key=lambda x:(x[5],x[0]),reverse=True) # 按照删除日期排序，删除日期相同时，按照项目号排序
    return delete_list 

def Generate_email_file(delete_list,outdir):  
    '''
    功能: 根据需要提醒的项目列表（已经按照日期和项目号排序），输出一个文件到输出路径，用于下一步发送邮件。outdir+'/email.txt'
    返回值：需要提醒的项目文件
    注意：输出的文件是按照项管整理的。
    '''
    persons = list(set([i[2] for i in delete_list]))
    out_file = outdir+'/email.txt'
    OUTfile = open(out_file,'w')
    title = '项目编号\t项目名称\t项目管理\t项目管理邮箱\t目录创建日期\t预计删除日期\t本地路径\n'
    OUTfile.write(title)
    write_list = [list(t) for t in set(tuple(_) for _ in delete_list)]
    write_list_sort = sorted(write_list,key=lambda x:(x[2],x[5],x[0]))  # 再次排序 去重之后就乱序了。,reverse=True

    for info in write_list_sort:
        content = "\t".join(info)+"\n"
        OUTfile.write(content)
    OUTfile.close()
    return out_file

def Send_email(config, email_ini, email_file, person, now_time):
    '''
    功能：发送邮件信息给每一个项管
    返回值：发送成功的标志文件 *.sign文件,发送邮件的命令
    注意事项: 失败和成功都会记录在日志文件里，成功之后会生成标志文件，失败则不会有。（后续衔接机器人，提醒发送失败的邮件。）
    '''
    PYTHON = config.get_value('PYTHON')
    email_prepare = Email(email_ini, email_file , person)
    email_dir, email_ini_out = email_prepare.email_info_write(now_time)
    send_qc_email_cmd = '{PYTHON} {bindir}/send_email.py -i {email_ini_out}'.format(PYTHON=PYTHON, bindir=bindir+'/../lib', email_ini_out=email_ini_out )
    run_state = myrun(send_qc_email_cmd)
    sign = email_file + ".sign"
    if run_state :
        my_log.info('数据超期删除提醒邮件  发送成功！！！')
        os.system("touch {}".format(sign))
    else:
        my_log.info('数据超期删除提醒邮件  发送失败，请检查原因！！！')
        my_log.info("可以尝试重新运行以下失败的命令： ".format( send_qc_email_cmd))
    return sign,send_qc_email_cmd

def main():
        parser=argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        epilog='author:\t{0}\nmail:\t{1}\n{2}'.format(__author__,__mail__,__doc__))
        parser.add_argument('-t','--test',help='test',dest='test',action='store_true')
        parser.add_argument('-e','--email',help='email Filter_Delete or not',dest='email',default=True)
        parser.add_argument('-s','--scan',help='scan release status',dest='scan',default=True)
        parser.add_argument('-o','--outdir',help='outdir',dest='outdir',required=True)
        args=parser.parse_args()
        config_dir = '{}/../config/'.format(bindir) 
        

        if args.test :
            table = 'tb_deletion_info_test'
            email_ini = '{}/pipe_email_test.ini'.format(config_dir)
            lims_config_file = config_dir+"/config_test.txt"
            config = Config( config_dir+"/config_test.txt" )  # 主要是 扫描的路径不同
        else:
            table = 'tb_deletion_info'
            email_ini = '{}/pipe_email.ini'.format(config_dir)
            lims_config_file = config_dir+"/config.txt"
            config = Config( config_dir+"/config.txt" )

        WarningDay = int(config.get_value('WarningDay'))  # 21，超期提醒时间，从过滤目录创建开始，经过多少天要提醒项管
        TimeFormat = config.get_value('TimeFormat')  # 时间格式，更新数据库时间时的格式
        DeleteDay = int(config.get_value('DeleteDay')) # 30，从过滤目录创建开始，经过多少天要删除该目录
        
        ReleaseDay = int(config.get_value('ReleaseDay')) # 15，硬盘交付超过多少天就会被判定为已经交付完成
        url = config.get_value('robot_url') # 机器人微信群的链接


        if args.scan == True :
            Scan_Status(lims_config_file,table,config,config_dir,TimeFormat,ReleaseDay)
            print("交付完成情况已更新！")

        if args.email == True :
            now_time = datetime.datetime.now().strftime("%Y-%m-%d") # 固定的，适配每天一封
            outdir = args.outdir+'/'+now_time  # 每天有一个目录，该目录下按照项管生成各自的提醒文件。
            mkdir(outdir)
            # 1、扫描本地路径，并存储为文件。如果已经存在，则不再扫描，对应发邮件失败后重跑。
            local_dir_file = outdir + "/all_local_dir_file.txt"
            if os.path.exists(local_dir_file) and os.path.getsize(local_dir_file):
                my_log.info("{}  文件已经存在且不为空，直接使用该文件筛选后续提醒项目 ".format(local_dir_file))
            else:
                CMD = config.get_value('Scan_CMD')
                Get_local_dir(CMD, local_dir_file)

            # 2、根据扫描的结果判断要提醒的项目列表
            delete_list = Scan_local_dir(WarningDay,lims_config_file,table,local_dir_file,TimeFormat,DeleteDay)
            

            # 3、根据项管名字和日期将列表排序，并生成要提醒的项目列表文件
            email_file = Generate_email_file(delete_list,outdir)

            # # 4、文件转码，发送项管群邮箱
            person = "Everyone"
            out_file = email_file.replace(".txt",'.{0}.xls'.format(now_time))
            CMD_convert = '/usr/bin/iconv -f utf-8 -t gbk {0} -o {1}'.format(email_file,out_file)

            if os.system(CMD_convert) == 0:
                # print(CMD_convert)
                print("文件转码成功，转码之后的文件是 {}".format(out_file))
            else:
                email_config = configparser.ConfigParser()
                email_config.read(email_ini)
                receiver = email_config['Filter_Delete']['Receiver']
                print("文件转码失败，请手动转码  {0} \n 并发送文件给以下邮箱， {1}".format(CMD_convert,receiver))
            if os.path.exists(out_file) and os.path.getsize(out_file):
                sign,cmd = Send_email( config, email_ini, out_file, person, now_time)
                n = 0
                while n < 3: # 如果失败，会再尝试3次
                    if os.path.exists(sign):
                        my_log.info("邮件已经发送成功，不再尝试")
                        break
                    else:
                        sign,cmd = Send_email( config, email_ini, out_file, person, now_time)
                        n += 1
                if n >2 :
                    my_log.info("邮件发送已经尝试4次， 皆失败，不再尝试 ，请核查原因 ")
                    Robots("【过滤项目数据超期删除提醒失败！！！请检查情况！】",out_file,cmd,url)

if __name__=="__main__":
        main()
