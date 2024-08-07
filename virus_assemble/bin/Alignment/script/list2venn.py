import sys
import re
import os
import argparse
from itertools import permutations, combinations
bindir = os.path.abspath(os.path.dirname(__file__))

__doc__='此脚本用于整理绘制韦恩图的前期数据，并生成对应out/InterGene下的文件\n使用示例：\npython3 list2venn.py -i example.xls -f "#436688,#7BBCFA,#6398CB,#77B4F1,#699ED2" -r "#000000,#000000,#000000,#000000,#000000" -a 0.6 -l solid -o ./output'
__author__='liaorui'
__mail__='ruiliao@genome.cn'
__date__='2021-04-16'

########获取每个文件的list
def get_list(input):
	name = []
	with open(input,'r') as f:
		for line in f:
			name.append(line.split("\t")[0].replace("\n",""))
	f.close()
	return name

#########输入为单个文件，将每一列转为一个txt文件
def get_file(input,num,output,title):
	file=open(input,'r')
	line=file.readlines()[1:]
	for j in range(len(title)):
		name = []
		with open(output+'/{0}.txt'.format(title[j]),'w') as out:
			for i in line:
				if i.split("\t")[j] == '':continue
				out.write("\t".join([str(i.split("\t")[j].replace('\n',''))])+"\n")
		out.close()
	file.close()

###########获取每一列的数目
def get_num(output,title):
	draw_count = {}
	for i in title:
		item = []
		with open(output+'/{0}.txt'.format(i),'r') as f:
			for line in f:
				if line.split('\t')[0] not in item:
					item.append(line.split('\t')[0])
		draw_count[i] = len(item)
	return draw_count

############写入excel文件中
def write_excel(list,output,type,name):
	with open(output+'/InterGene/{0}/{1}.xls'.format(type,name),'w') as out:
		out.write('name\n')
		for i in list:
			i.replace("\n","")
		out.write('\n'.join(list))
	out.close()

###########根据排列组合获取每一个组合，对每一个组合获取共同的部分
def get_group(title,output,type):
	i = 2
	draw_count = {}
	while i <= len(title):
		com = list(combinations(title,i))
		for group in range(len(com)):
			sample = 1
			common = get_list(output+'/{0}.txt'.format(com[group][0]))
			while sample < len(com[group]):
				file = get_list(output+'/{0}.txt'.format(com[group][sample]))
				common = list(set(file) & set(common))
				sample += 1
			name = 'VS'.join(com[group])
			out2 = write_excel(common,output,type['{0}'.format(len(com[group]))],name)
			draw_count[name] = len(set(common))
		i += 1
	return draw_count

#########对输入的color颜色进行处理，如果输入多了则截取，如果少了则进行循环获取需要的颜色数目
def color_num(color,title):
	color_num = color.split(",")
	if len(color_num) == len(title):
		color_new = ','.join(color_num)
	elif len(color_num) > len(title):
		color_new = ','.join(color_num[:len(title)])
	elif len(color_num) < len(title):
		color_plus = color_num*5
		color_new = ','.join(color_num[:len(title)])
	return color_new

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file,split with ","',dest='input')
	parser.add_argument('-f','--fill_color',help='fill color,split with ","',dest='fill',default='#004DA1,#F7CA18,#4ECDC4,#F9690E,#B35AA5')
	parser.add_argument('-r','--rim_color',help='rim color,split with ","',dest='rim',default='NA')
	parser.add_argument('-a','--alpha',help='transparency',dest='alpha',default='6')
	parser.add_argument('-l','--lty',help='lty,please choose one from solid,dotted,blank',dest='lty',default='solid')
	parser.add_argument('-o','--outdir',help='outdir',dest='outdir')
	args=parser.parse_args()

####获取文件数目以及title
	type = {'2':'pairwise','3':'triple','4':'quad','5':'quintuple'}
	outdir = os.path.abspath(args.outdir)
	if os.path.isdir(outdir+'/InterGene/') == True:pass
	else:
		os.makedirs(outdir+'/InterGene/')

#########获取title
	file_num = len(args.input.split(","))
	if file_num == 1:
		title = open(args.input).readline().replace('\n','').split('\t')
		if len(set(title)) != len(title):
			print("有重复的名称，烦请检查")
			sys.exit()
		else:
			out1 = get_file(args.input,len(title),outdir,title)
	else:
		title = []
		for i in range(len(args.input.split(","))):
			title_name = os.path.splitext(os.path.basename(args.input.split(",")[i]))[0]
			if title_name in title:
				print("有重复的名称，烦请检查")
				sys.exit()
			else:
				title.append(title_name)
				os.system('cat {0}|cut -f1 >{1}/{2}.txt'.format(args.input.split(",")[i],outdir,title_name))
	print("获取样本名称为：{0}".format(','.join(title)))

#########建立目录
	file_type = type['{0}'.format(len(title))]
	i = 2
	while i <= len(title):
		if os.path.isdir(outdir+'/InterGene/'+type['{0}'.format(i)]) == True:pass
		else:
			os.makedirs(outdir+'/InterGene/'+type['{0}'.format(i)])
		i += 1

#####进行交并集处理，并输出至InterGene文件夹
	if len(title) <= 1 or len(title) > 5:
		print("Error：韦恩图需要2个及以上类型，或者请检查文件格式：每一列为一个类型")
		sys.exit()
	else:
		com = get_group(title,outdir,type)
		uniq = get_num(outdir,title)
		with open(outdir+'/draw_venn.xls','w') as out2:
			for i in range(len(title)):
				out2.write('\t'.join([title[i],str(uniq[title[i]])])+'\n')
#				out2.write('\t'.join([title[i],str(uniq[title[i]]),str(len(title[i]))])+'\n')
			j = 2
			while j <= len(title):
				combine = list(combinations(title,j))
				for group in range(len(combine)):
					com_name = 'VS'.join(combine[group])
					out2.write('\t'.join([com_name,str(com[com_name])])+'\n')
				j += 1
#			for key in com:
#				out2.write('\t'.join([key,str(com[key]),str(len(key))])+'\n')
		out2.close()


#####开始绘图
#	os.system('cat {0}/draw_venn.tmp.xls|sort -k3 -n|cut -f 1-2 >{0}/draw_venn.xls'.format(outdir))
	for i in range(len(args.fill.split(","))):
		if args.fill.split(",")[i].startswith('#') != True or len(args.fill.split(",")[i]) != 7:
			print("Warning：颜色并非以#开头的十六进制颜色代码，可能无法绘图")
	for i in range(len(args.rim.split(","))):
		if args.rim.split(",")[i].startswith('#') != True or len(args.rim.split(",")[i]) != 7:
			print("Warning：颜色并非以#开头的十六进制颜色代码，可能无法绘图")
	fill_color = color_num(args.fill,title)
	rim_color = color_num(args.rim,title)
	print('执行脚本：Rscript {0}/venn_r/draw.{1}.venn.r {2}/draw_venn.xls {2}/venn.pdf "{3}" "{4}" {5} {6}'.format(bindir,type['{0}'.format(len(title))],outdir,fill_color,rim_color,args.alpha,args.lty))
	os.system('Rscript {0}/venn_r/draw.{1}.venn.r {2}/draw_venn.xls {2}/venn.pdf "{3}" "{4}" {5} {6}'.format(bindir,type['{0}'.format(len(title))],outdir,fill_color,rim_color,args.alpha,args.lty))
	os.system('/usr/bin/convert {0}/venn.pdf {0}/venn.png'.format(outdir))
	os.system('rm {0}/*.txt'.format(outdir))
#	os.system('rm {0}/*.xls'.format(outdir))

if __name__ == '__main__':
	main()
