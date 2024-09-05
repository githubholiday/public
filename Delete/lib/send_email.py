'''
程序说明
'''
import os
import sys
import argparse
bindir = os.path.abspath(os.path.dirname(__file__))
sys.path.append( bindir )
import myemail

__author__ = 'suyanxun'
__mail__ = 'yanxunsu@genome.cn'

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:	{0}	mail:	{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--infile',help='infile',dest='infile',nargs = '+',required=True)
	args=parser.parse_args()
	filelist = []
	filelist.extend( args.infile )
	email = myemail.Email( filelist )
	email.send_email()

if __name__=="__main__":
	main()
