import threading
import time
import queue
import subprocess
import sys
import os
import re
import signal

def myTime():
	return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) 


class MyThread(threading.Thread):
	def __init__(self ,  thread_name , logfile, cmd , major , depend = [] ):
		threading.Thread.__init__(self)
		#self.lock = lock
		self.name = thread_name
		self.cmd = cmd
		self.logfile = logfile
		self.major = major
		self.status = 'waiting' 
		self.depend = depend
		self.bool_depend_job_finish = True
	def run(self):
		#self.lock.acquire()
		self.logfile.write('[start]: {0} start at {1}\n'.format(self.name , myTime()))
		self.logfile.flush()
		print(self.cmd)
		if subprocess.call(self.cmd , shell=True) != 0 : 
			self.logfile.write('[break]: {0} break down at {1}\n'.format(self.name , myTime()))
			self.logfile.flush()
			self.status = 'break'
			#sys.exit()
		#if self.status == True:
		else:
			self.logfile.write('[finish]: {0} finish at {1}\n'.format(self.name , myTime()))
			self.logfile.flush()
			self.status = 'finish'
		#self.lock.release()

def run(jobs , logfile):
	unsubmit_depend_job = []
	unfinish_job = []
	for order in sorted(jobs):
		break_list = []
		m = []
		for a_job in jobs[order] + unsubmit_depend_job :
			if a_job.status == 'finish' : continue
			bool_finish = True
			for pre_job in a_job.depend:
				if not pre_job.status == 'finish': 
					if pre_job.status == 'waiting':
						#time.sleep(300)
						if not a_job in unsubmit_depend_job:
							unsubmit_depend_job.append(a_job)
						bool_finish = False
					elif pre_job.status == 'break' :
						pre_job.logfile.write('#{0} break, Sorry,goodbye\n'.format(pre_job.name))
						sys.exit("{0} break, Sorry,goodbye".format(pre_job.name))
						#sys.exit("Sorry,goodbye")
			a_job.bool_depend_job_finish = bool_finish
			if a_job.bool_depend_job_finish : 
				a_job.start()  
				unfinish_job.append(a_job)
			if a_job.major == True : m.append(a_job)
		for a_job in m:
			a_job.join()
			time.sleep(60)
		for a_job in jobs[order]:
			if a_job == '' : continue
			if a_job.status == 'break' :
				a_job.logfile.write('#{0} break, Sorry,goodbye\n'.format(a_job.name))
				if a_job in m :
					break_list.append(a_job.name)
				else :
					pass
		if len(break_list) > 0 : 
			sys.exit("{0} break, Sorry,goodbye\n".format("|".join(break_list)))

	while unsubmit_depend_job :
		a_job = unsubmit_depend_job.pop(0)
		#a_job.logfile.write("unsubmit:{0}\n".format("\t".join( [i.name for i in unsubmit_depend_job])))
		if a_job.status == 'finish' : continue
		bool_finish = True
		for pre_job in a_job.depend:
			if not pre_job.status == 'finish': 
				if pre_job.status == 'waiting':
						#time.sleep(300)
					if not a_job in unsubmit_depend_job:
						unsubmit_depend_job.append(a_job)
					bool_finish = False
				elif pre_job.status == 'break' :
					pre_job.logfile.write('#{0} break, Sorry,goodbye\n'.format(pre_job.name))
					sys.exit("{0} break, Sorry,goodbye".format(pre_job.name))
						#sys.exit("Sorry,goodbye")
		a_job.bool_depend_job_finish = bool_finish
		if a_job.bool_depend_job_finish : 
			a_job.start()  
			unfinish_job.append(a_job)
		time.sleep(30)

	while unfinish_job:
		a_job = unfinish_job.pop(0)
		#a_job.logfile.write("unfinish:{0}\n".format("\t".join([i.name for i in unfinish_job])))
		if a_job == '' : continue
		if a_job.status == 'finish':
			continue
		else:
			unfinish_job.append(a_job)
		if a_job.status == 'break' :
			a_job.logfile.write('#{0} break, Sorry,goodbye\n'.format(a_job.name))
		time.sleep(30)

def RemoveFinish(jobs, finish):
	for order in jobs:
		for count, a_job in enumerate(jobs[order]):
			if a_job.name in finish:
				a_job.status='finish'
	return jobs

def ReadLog(logfile):
	pat1 = re.compile('^\s+$')
	finish = []
	if not os.path.isfile(logfile):
		return finish
	else:
		with open(logfile,'r') as f_file:
			for line in f_file:
				if line.startswith('#') or re.search(pat1,line):continue
				tmp=line.rstrip().split()
				if line.startswith('[finish]:'):
					finish.append(tmp[1])
		return finish

def signal_term_handler(signal , frame , logfile ):
	logfile.write('[break]: system was break down at {0}\n'.format( myTime()))
	sys.exit(0)




