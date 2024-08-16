
'''
Parameters:
	-i , --input : 输入的项目配置文件
	-b , --bin : 程序调用的相关程序路径，与流程配置文件中的${BIN}对应
	-t , --thread: 线程数，如果定义了这个，那么会覆盖流程配置文件中的CPU
	-q , --queue : 指定的队列，如果定义了这个，则会覆盖流程配置文件中的队列
	-P , --P : 指定的qsub的Project队列，默认是none，即为普通项目，如果需要加急，可以设置为vip（全小写）
	-o , --outdir: 输出目录，与流程配置文件中的$(OUTDIR)对应
	-name,--name : 任务名称
	-j,   --jobid: 任务前缀，默认为name
	-r,   --run  : 是否自动投递任务，默认为不投递任务，但会在${OUTDIR}/shell中生成脚本，可以进行检查
	-c,   --continue: 在qsub下有效（设置了-r，但不设置-n），如果某一步分析中没有完成全部任务，如果不指定则从头运行该步所有任务，指定则完成该步剩余未完成任务
	-a,   --add   : 是否只运行加测的样品，如果指定则只运行加测的样品。注意的是，-a不可与-c同时使用。如果使用-a，需要将之前的log文件删除干净，否则不会运行任务；
	-quota, --quota : 分析目录的配额，是之前找文明申请的大小。默认是1000G，请根据实际情况进行调整。
说明：
	Q: 支线任务断了咋办？
	A: 支线任务断了，主线任务会继续运行，而整个流程不受影响。但是查看show_process的时候，可以看到break的状态。而break的等级是高于run和hold的，因此，需要自行判定程序是否完成。

	Q: 主线任务断了咋办？
	A: 主线断了，那么就断了，需要重新投递。
	   如果主线任务断了，而支线任务没有完成，那么主程序会等候支线任务完成才会实现最终的退出，所以会导致程序一直在运行的假象。这个时候有两个处理办法：
	   1. 把所有进程都杀掉，然后重新投递
	   2. 只杀掉主进程（pipeline.py)，并且在log.txt中人为添加支线任务的finish标识（防止再次运行pipeline.py时重新投递），之后重新运行该任务。
	Q:如何识别加测样品？
	A:当某一块出现两次或者多次，那么最后一次出现的内容作为加测项

	Q:如何监控项目运行状态？
	A:运行/annoroad/bioinfo/PMO/liutao/pipeline_generate/bin/current/show_process.py会显示项目的状态。项目运行状态分为running, break, plan, end 四种。其中running表示正在运行，break表示中断，plan表示准备运行，end表示运行完成,hold表示磁盘不够，任务挂起。

	Q: 发现任务状态是break，该咋办？
	A: 当发现任务状态是break的时候，首先需要确定break掉的任务是否是主线任务。
	   如果是主线任务，如果主程序自然退出(在top或者ps的时候没有发现pipeline.py），则可以重新投递任务；
	                   如果没有自然退出（主程序pipeline.py还在运行），那么可能是有之前的支线任务未完成，可以参照Q2来进行操作；
	   如果是支线任务，如果查看log.txt发现主线任务全部完成或者正常运行，如果时间允许，可以等待所有任务完成后再投递任务；
	                                                                    如果加急，可以把该步对应的sh文件修改后，手动投递该任务；
	                   如果主线任务也断掉了，那么修改脚本后，重新投递所有任务。

	Q:监控项目的记录文件在哪？
	A:程序会在您的home目录下，生成一个记录文件，路径为~/.mission/.pipeline.log，记录了每个项目的分析目录。如果不想显示某个项目，可以对相应的行进行删除或者编辑。

	Q:如果程序断了，咋办？
	A:如果程序由于各种因素中断了，仔细检查脚本，如果脚本没错，确定只是中断，那么重新运行一次之前的脚本，默认会把断掉的模块全部重头运行；如果不想将该模块内已经完成的样品重新运行，可以加上-c参数，那么会只运行没有运行成功的样品。

	Q:如何发现配额不够？
	A:如果配额不够了，会把所有的任务挂起，使用show_process查看时，会发现Hold状态；或者使用qstat的时候，会发现hqw，hr，ht等，或者没有任务在运行（因为配额不够，程序会自动不投递任务）

	Q: 配额不够了，咋办？
	A: 第一，找文明修改配额 ； 第二，修改相应的sh.*.log文件，加入一行DISK_QUOTA	**G ;之后，程序会自动的release。但需要注意的是，因为之前在程序里设置了较小的配额，所以之后每一步都会被hold。所以需要之后每个log文件都加上 DISK_QUOTA	**G，来每次进行更改；或者杀掉重新来。 
	   或者删除文件来释放空间，这样的话，后面可以不用修改就可以运行。

	Q: 如何杀掉程序？
	A: 1. 杀掉所有的子进程 守护进程qsub_sge.pl，否则的杀掉的任务会重新投递；
	   2. 杀掉所有的任务 qdel掉

	Q: 如何精准的杀掉守护进程？
	A: 1. 在重新投递之前，使用ps -f -u name |cat 然后仔细的判别，获得进程ID
	   2. 查看shell后面的数字，在sh 和log直接的数字，是进程ID，使用kill -9可以杀掉
更新说明：
	2015-8-17
	1.之前一个target断掉后，需要对这个target里的所有任务重新投递；目前加上了-c 参数，可以选择跳过已运行完毕的任务
	2015-9-1
	1. 加入quota参数，如果目录配额不够，会自动挂起任务；
	2. 加入sh.*.log文件中的节点的记录
	3. 加入了maxcycle，可以在pipeline.py对maxcycle中进行修改
	4. 修改支线任务断掉，不影响主线任务整体运行。
	2015-10-8
	1. 添加了每个任务投递时的线程数，在最初的config文件中使用Thread控制，默认为1.
	2016-03-26
	1. 在shell下面的log文件中，添加了使用程序版本的记录，方便以后升级时，可以进行查找。
'''
#! /usr/bin/env python3
# -*- coding: utf-8 -*-  
import argparse
import sys
import os
import re
bindir = os.path.abspath(os.path.dirname(__file__))
sys.path.append('{0}/lib'.format(bindir))
import parseConfig
import JobGuard

__author__='Liu Tao'
__mail__= 'taoliu@annoroad.com'

pat1=re.compile('^\s+$')

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',type=open,required=True)
	parser.add_argument('-b','--bin',help='bin dir',dest='bin',default=os.path.dirname(bindir))
	parser.add_argument('-t','--thread',help='thread number ',dest='thread',type=int)
	parser.add_argument('-q','--queue',help='computer queue name ',dest='queue')
	parser.add_argument('-P','--P',help='Project in qsub ',dest='qsub_P' , default='none')
	parser.add_argument('-o','--outdir',help='output file',dest='outdir',required=True)
	parser.add_argument('-indir','--indir',help='indir data file',dest='indir',default='')
	parser.add_argument('-name','--name',help='project name',dest='name',required=True)
	parser.add_argument('-j','--jobid',help='job id prefix',dest='jobid',default='')
	parser.add_argument('-r','--run',help='run script file',dest='run',action='store_true')
	parser.add_argument('-nc','--noncontinue',help='continue unfinish job in each shell',dest='noncontinues',action='store_true')
	#parser.add_argument('-n','--nohup',help='qsub or nohup mission',dest='nohup',action='store_true')
	parser.add_argument('-quota','--quota',help='disk quota ',dest='quota',default = '1000G')
	parser.add_argument('-a','--add',help='add sequencing sample process, True -- only run added sequence sample ,False* -- run all samples',dest='add',action='store_true')
	args=parser.parse_args()

	OUTDIR = parseConfig.getab(args.outdir)
	INDIR = parseConfig.getab(args.indir)
	BIN=os.path.realpath(args.bin)
	LOGFILE = '{0}/log.txt'.format(OUTDIR)
	
	job_not_continue = '  '
	if args.noncontinues : job_not_continue = '  -nc  '

	if args.jobid == '' : args.jobid = args.name
	config ,para, db , orders  = parseConfig.ReadConfig(args.input)
	shell_dir = '{0}/shell/'.format(OUTDIR)
	parseConfig.makedir(shell_dir)
	logfile = '{0}/log.txt'.format(shell_dir)
	finish_obj = JobGuard.ReadLog(logfile)
	log = open(logfile,'a')
	log.write('#pipeline version : {0}\n'.format(BIN))
	guard_script = '{0}/guard.py'.format(shell_dir)
	job_list = {}
	all_job_list = {}

	run_sample = config['sample']
	pre_job_count = 0 
	if args.add :
		run_sample , pre_job_count = parseConfig.chooseSamples(run_sample, orders['sample'])
	#print(run_sample)
	if len(run_sample) > 0 :
		shsh = '{0}/1_0_MarkerAnno.sh'.format(shell_dir)
		with open(shsh, 'w') as f_out:
			cmds = []
			cpu = parseConfig.cpu(args.thread , len(run_sample), 'N' )
			queue = parseConfig.queue(args.queue , 'sci.q')
			for sample in sorted(run_sample):
				cmds.append('make -f {BIN}/bin/anno.mk log_file={LOGFILE} sample={sample[0]} infile={sample[1]} outdir={sample[2]} relation_file={para[Para_relation_file]} anno_file={para[Para_anno_file]} Anno'.format(para=para , sample =sample ,OUTDIR=OUTDIR, BIN=BIN,db=db,LOGFILE=LOGFILE, INDIR=INDIR) )
			f_out.write("\n".join(sorted(set(cmds))) + '\n')
		if not False :
			if "" == '':
				if "" == '':
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 {0}'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 --mount / --sif   {0}'.format(shsh , bindir, cpu)
			else:
				if "" == '':
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 {0}"'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 --mount / --sif   {0}"'.format(shsh , bindir, cpu)
		else:
			if "" == '':
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}MarkerAnno -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}MarkerAnno -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
			else :
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}MarkerAnno -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}MarkerAnno -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
		if len( [] ) == 0 : 
			a_thread = JobGuard.MyThread('1_0_MarkerAnno' , log , a_cmd, False)
		else : 
			a_thread = JobGuard.MyThread('1_0_MarkerAnno' , log , a_cmd, False , [ all_job_list[ i ] for i in [] ])
		all_job_list[ 'MarkerAnno' ] = a_thread
		if not int(1) in job_list: job_list[int(1)] = []
		job_list[int(1)].append(a_thread)
	else:
		print("{0} is empty".format("config['sample']"))

	run_sample = config['sample']
	pre_job_count = 0 
	if args.add :
		run_sample , pre_job_count = parseConfig.chooseSamples(run_sample, orders['sample'])
	#print(run_sample)
	if len(run_sample) > 0 :
		shsh = '{0}/2_0_UpGO.sh'.format(shell_dir)
		with open(shsh, 'w') as f_out:
			cmds = []
			cpu = parseConfig.cpu(args.thread , len(run_sample), 'N' )
			queue = parseConfig.queue(args.queue , 'sci.q')
			for sample in sorted(run_sample):
				cmds.append('make -f {BIN}/bin/anno.mk log_file={LOGFILE} outdir={sample[2]}/up/ up_or_down=Up infile={sample[2]}/$(sample).csv GetList\nmake -f {BIN}/bin/anno.mk log_file={LOGFILE} gene_list={sample[2]}/up/Up.gene.xls go_dir={sample[2]}/up/GO sample={sample[0]} go={para[Para_go_list]} GO'.format(para=para , sample =sample ,OUTDIR=OUTDIR, BIN=BIN,db=db,LOGFILE=LOGFILE, INDIR=INDIR) )
			f_out.write("\n".join(sorted(set(cmds))) + '\n')
		if not False :
			if "" == '':
				if "" == '':
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 {0}'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 --mount / --sif   {0}'.format(shsh , bindir, cpu)
			else:
				if "" == '':
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 {0}"'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 --mount / --sif   {0}"'.format(shsh , bindir, cpu)
		else:
			if "" == '':
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}UpGO -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}UpGO -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
			else :
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}UpGO -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}UpGO -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
		if len( ['MarkerAnno'] ) == 0 : 
			a_thread = JobGuard.MyThread('2_0_UpGO' , log , a_cmd, False)
		else : 
			a_thread = JobGuard.MyThread('2_0_UpGO' , log , a_cmd, False , [ all_job_list[ i ] for i in ['MarkerAnno'] ])
		all_job_list[ 'UpGO' ] = a_thread
		if not int(2) in job_list: job_list[int(2)] = []
		job_list[int(2)].append(a_thread)
	else:
		print("{0} is empty".format("config['sample']"))

	run_sample = config['sample']
	pre_job_count = 0 
	if args.add :
		run_sample , pre_job_count = parseConfig.chooseSamples(run_sample, orders['sample'])
	#print(run_sample)
	if len(run_sample) > 0 :
		shsh = '{0}/2_1_DownGO.sh'.format(shell_dir)
		with open(shsh, 'w') as f_out:
			cmds = []
			cpu = parseConfig.cpu(args.thread , len(run_sample), 'N' )
			queue = parseConfig.queue(args.queue , 'sci.q')
			for sample in sorted(run_sample):
				cmds.append('make -f {BIN}/bin/anno.mk log_file={LOGFILE} outdir={sample[2]}/down/ up_or_down=Down infile={sample[2]}/$(sample).csv GetList\nmake -f {BIN}/bin/anno.mk log_file={LOGFILE} gene_list={sample[2]}/down/Down.gene.xls go_dir={sample[2]}/down/GO sample={sample[0]} go={para[Para_go_list]} GO'.format(para=para , sample =sample ,OUTDIR=OUTDIR, BIN=BIN,db=db,LOGFILE=LOGFILE, INDIR=INDIR) )
			f_out.write("\n".join(sorted(set(cmds))) + '\n')
		if not False :
			if "" == '':
				if "" == '':
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 {0}'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 --mount / --sif   {0}'.format(shsh , bindir, cpu)
			else:
				if "" == '':
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 {0}"'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 --mount / --sif   {0}"'.format(shsh , bindir, cpu)
		else:
			if "" == '':
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}DownGO -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}DownGO -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
			else :
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}DownGO -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}DownGO -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
		if len( ['MarkerAnno'] ) == 0 : 
			a_thread = JobGuard.MyThread('2_1_DownGO' , log , a_cmd, False)
		else : 
			a_thread = JobGuard.MyThread('2_1_DownGO' , log , a_cmd, False , [ all_job_list[ i ] for i in ['MarkerAnno'] ])
		all_job_list[ 'DownGO' ] = a_thread
		if not int(2) in job_list: job_list[int(2)] = []
		job_list[int(2)].append(a_thread)
	else:
		print("{0} is empty".format("config['sample']"))

	run_sample = config['sample']
	pre_job_count = 0 
	if args.add :
		run_sample , pre_job_count = parseConfig.chooseSamples(run_sample, orders['sample'])
	#print(run_sample)
	if len(run_sample) > 0 :
		shsh = '{0}/2_2_AllGO.sh'.format(shell_dir)
		with open(shsh, 'w') as f_out:
			cmds = []
			cpu = parseConfig.cpu(args.thread , len(run_sample), 'N' )
			queue = parseConfig.queue(args.queue , 'sci.q')
			for sample in sorted(run_sample):
				cmds.append('make -f {BIN}/bin/anno.mk log_file={LOGFILE} gene_list={sample[2]}/$(sample).csv go_dir={sample[2]}/all/GO sample={sample[0]} go={para[Para_go_list]} GO'.format(para=para , sample =sample ,OUTDIR=OUTDIR, BIN=BIN,db=db,LOGFILE=LOGFILE, INDIR=INDIR) )
			f_out.write("\n".join(sorted(set(cmds))) + '\n')
		if not False :
			if "" == '':
				if "" == '':
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 {0}'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 --mount / --sif   {0}'.format(shsh , bindir, cpu)
			else:
				if "" == '':
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 {0}"'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 --mount / --sif   {0}"'.format(shsh , bindir, cpu)
		else:
			if "" == '':
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}AllGO -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}AllGO -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
			else :
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}AllGO -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}AllGO -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
		if len( ['MarkerAnno'] ) == 0 : 
			a_thread = JobGuard.MyThread('2_2_AllGO' , log , a_cmd, False)
		else : 
			a_thread = JobGuard.MyThread('2_2_AllGO' , log , a_cmd, False , [ all_job_list[ i ] for i in ['MarkerAnno'] ])
		all_job_list[ 'AllGO' ] = a_thread
		if not int(2) in job_list: job_list[int(2)] = []
		job_list[int(2)].append(a_thread)
	else:
		print("{0} is empty".format("config['sample']"))

	run_sample = config['sample']
	pre_job_count = 0 
	if args.add :
		run_sample , pre_job_count = parseConfig.chooseSamples(run_sample, orders['sample'])
	#print(run_sample)
	if len(run_sample) > 0 :
		shsh = '{0}/2_3_UpKEGG.sh'.format(shell_dir)
		with open(shsh, 'w') as f_out:
			cmds = []
			cpu = parseConfig.cpu(args.thread , len(run_sample), 'N' )
			queue = parseConfig.queue(args.queue , 'sci.q')
			for sample in sorted(run_sample):
				cmds.append('make -f {BIN}/bin/anno.mk log_file={LOGFILE} outdir={sample[2]}/up/ up_or_down=Up infile={sample[2]}/$(sample).csv GetList\nmake -f {BIN}/bin/anno.mk log_file={LOGFILE} gene_list={sample[2]}/up/Up.gene.xls go_dir={sample[2]}/up/GO sample={sample[0]} kegg={para[Para_kegg_list]} category=animal KEGG'.format(para=para , sample =sample ,OUTDIR=OUTDIR, BIN=BIN,db=db,LOGFILE=LOGFILE, INDIR=INDIR) )
			f_out.write("\n".join(sorted(set(cmds))) + '\n')
		if not False :
			if "" == '':
				if "" == '':
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 {0}'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 --mount / --sif   {0}'.format(shsh , bindir, cpu)
			else:
				if "" == '':
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 {0}"'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 --mount / --sif   {0}"'.format(shsh , bindir, cpu)
		else:
			if "" == '':
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}UpKEGG -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}UpKEGG -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
			else :
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}UpKEGG -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}UpKEGG -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
		if len( ['MarkerAnno'] ) == 0 : 
			a_thread = JobGuard.MyThread('2_3_UpKEGG' , log , a_cmd, False)
		else : 
			a_thread = JobGuard.MyThread('2_3_UpKEGG' , log , a_cmd, False , [ all_job_list[ i ] for i in ['MarkerAnno'] ])
		all_job_list[ 'UpKEGG' ] = a_thread
		if not int(2) in job_list: job_list[int(2)] = []
		job_list[int(2)].append(a_thread)
	else:
		print("{0} is empty".format("config['sample']"))

	run_sample = config['sample']
	pre_job_count = 0 
	if args.add :
		run_sample , pre_job_count = parseConfig.chooseSamples(run_sample, orders['sample'])
	#print(run_sample)
	if len(run_sample) > 0 :
		shsh = '{0}/2_4_DownKEGG.sh'.format(shell_dir)
		with open(shsh, 'w') as f_out:
			cmds = []
			cpu = parseConfig.cpu(args.thread , len(run_sample), 'N' )
			queue = parseConfig.queue(args.queue , 'sci.q')
			for sample in sorted(run_sample):
				cmds.append('make -f {BIN}/bin/anno.mk log_file={LOGFILE} outdir={sample[2]}/down/ up_or_down=Down infile={sample[2]}/$(sample).csv GetList\nmake -f {BIN}/bin/anno.mk log_file={LOGFILE} gene_list={sample[2]}/down/Down.gene.xls go_dir={sample[2]}/down/GO sample={sample[0]} kegg={para[Para_kegg_list]} category=animal KEGG'.format(para=para , sample =sample ,OUTDIR=OUTDIR, BIN=BIN,db=db,LOGFILE=LOGFILE, INDIR=INDIR) )
			f_out.write("\n".join(sorted(set(cmds))) + '\n')
		if not False :
			if "" == '':
				if "" == '':
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 {0}'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 --mount / --sif   {0}'.format(shsh , bindir, cpu)
			else:
				if "" == '':
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 {0}"'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 2 --mount / --sif   {0}"'.format(shsh , bindir, cpu)
		else:
			if "" == '':
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}DownKEGG -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}DownKEGG -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
			else :
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}DownKEGG -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 2 --maxcycle 5 --quota {6}  --jobprefix {3}DownKEGG -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
		if len( ['MarkerAnno'] ) == 0 : 
			a_thread = JobGuard.MyThread('2_4_DownKEGG' , log , a_cmd, False)
		else : 
			a_thread = JobGuard.MyThread('2_4_DownKEGG' , log , a_cmd, False , [ all_job_list[ i ] for i in ['MarkerAnno'] ])
		all_job_list[ 'DownKEGG' ] = a_thread
		if not int(2) in job_list: job_list[int(2)] = []
		job_list[int(2)].append(a_thread)
	else:
		print("{0} is empty".format("config['sample']"))

	run_sample = config['sample']
	pre_job_count = 0 
	if args.add :
		run_sample , pre_job_count = parseConfig.chooseSamples(run_sample, orders['sample'])
	#print(run_sample)
	if len(run_sample) > 0 :
		shsh = '{0}/2_5_AllKEGG.sh'.format(shell_dir)
		with open(shsh, 'w') as f_out:
			cmds = []
			cpu = parseConfig.cpu(args.thread , len(run_sample), 'N' )
			queue = parseConfig.queue(args.queue , 'sci.q')
			for sample in sorted(run_sample):
				cmds.append('make -f {BIN}/bin/anno.mk log_file={LOGFILE} gene_list={sample[2]}/$(sample).csv go_dir={sample[2]}/all/GO sample={sample[0]} kegg={para[Para_kegg_list]} category=animal KEGG'.format(para=para , sample =sample ,OUTDIR=OUTDIR, BIN=BIN,db=db,LOGFILE=LOGFILE, INDIR=INDIR) )
			f_out.write("\n".join(sorted(set(cmds))) + '\n')
		if not False :
			if "" == '':
				if "" == '':
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 {0}'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 --mount / --sif   {0}'.format(shsh , bindir, cpu)
			else:
				if "" == '':
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 {0}"'.format(shsh , bindir, cpu)
				else :
					a_cmd = 'ssh  2> /dev/null "perl {1}/src/slurm/multi-process.pl -cpu {2} --lines 1 --mount / --sif   {0}"'.format(shsh , bindir, cpu)
		else:
			if "" == '':
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}AllKEGG -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}AllKEGG -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
			else :
				if "" == '':
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}AllKEGG -P {8}  {5} --queue {4} {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
				else :
					a_cmd = '{1}/src/slurm/slurm_sge --resource "p=1,vf=5G,h=" --maxjob {2} --lines 1 --maxcycle 5 --quota {6}  --jobprefix {3}AllKEGG -P {8}  {5} --queue {4}  -mount / -s  {0} --jobidstart {7} '.format(shsh ,bindir, cpu, args.jobid , queue , job_not_continue , args.quota , pre_job_count , args.qsub_P)
		if len( ['MarkerAnno'] ) == 0 : 
			a_thread = JobGuard.MyThread('2_5_AllKEGG' , log , a_cmd, False)
		else : 
			a_thread = JobGuard.MyThread('2_5_AllKEGG' , log , a_cmd, False , [ all_job_list[ i ] for i in ['MarkerAnno'] ])
		all_job_list[ 'AllKEGG' ] = a_thread
		if not int(2) in job_list: job_list[int(2)] = []
		job_list[int(2)].append(a_thread)
	else:
		print("{0} is empty".format("config['sample']"))

	home_dir = os.environ['HOME']
	parseConfig.makedir('{0}/.mission'.format(home_dir))
	if not os.path.isfile('{0}/.mission/.pipeline.log'.format(home_dir)):
		os.system('touch {0}/.mission/.pipeline.log'.format(home_dir))
	
	tag = 0 
	with open('{0}/.mission/.pipeline.log'.format(home_dir),'r') as super_log:
		for line in super_log:
			if line.startswith('#') or re.search(pat1,line):continue
			tmp = line.rstrip().split()
			if tmp[0] == args.name:
				tag = 1
				if tmp[1] == os.path.abspath(args.outdir):
					tag = 2
	
	if tag == 0 : 
		with open('{0}/.mission/.pipeline.log'.format(home_dir),'a') as super_log:
			super_log.write('{0}\t{1}\n'.format(args.name , os.path.abspath(args.outdir)))
	elif tag == 2 :
		print("[1;31;40m" + "Warings: {0} was existed already in your log file, please check it".format(args.name) + "[0m")
	elif tag == 1 :
		print("[1;31;40m" + "Warings: {0} was existed already in your log file,  and have different analysis directory , we should add this new dir at the end of log file ,please check it".format(args.name) + "[0m")
		with open('{0}/.mission/.pipeline.log'.format(home_dir),'a') as super_log:
			super_log.write('{0}\t{1}\n'.format(args.name , os.path.abspath(args.outdir)))
	
	print(JobGuard.dump_jobs(job_list))
	job_list = JobGuard.RemoveFinish(job_list,finish_obj)
	if args.run == True:
		JobGuard.run(job_list , log)
	log.close()

if __name__ == '__main__':
	main()
