import os
import sys
import argparse
import glob
import comm

__author__ = 'tuchengfang'
__mail__ = 'chengfangtu@genome.cn'

class MD5( object ) :
    def __init__(self, release_dir ) :
        self.release_dir = release_dir
        self.md5file = self.release_dir+'/md5.txt'
        #self.md5_repeat = self.check_repeat_md5()

    def check_repeat_md5( self ) :
        md5_dic = {}
        self.md5_line_count = 0
        with open( self.md5file, 'r') as md5txt :
            for line in md5txt :
                self.md5_line_count += 1
                tmp = line.rstrip().split('  ')
                md5_code = tmp[0]
                md5_file = tmp[1]
                if md5_code not in md5_dic :
                    md5_dic[md5_code] = []
                md5_dic[md5_code].append( md5_file )
        md5_repeat = ''
        for md5 in md5_dic :
            md5_f = md5_dic[md5]
            if len(md5_f) > 1 :
                for file in md5_f :
                    md5_repeat += '{0}  {1}\n'.format(md5, file)
        if md5_repeat :
            comm.color_print("md5有重复:{0},程序退出".format(md5_repeat))
            sys.exit(1)
        else:
            return md5_dic
        
    def compare_md5_sample( self, sample_list, md5_dic ):
        '''count the file number in release dir
        sample_list:以逗号连接的字符串
        md5_dic[md5_value] = [Cleandata/sample/sample_R1.fq.gz]
        '''
        md5_sample_list = []
        for md5_value in md5_dic:
            filename = md5_dic[md5_value][0]
            if '/' not in filename: continue
            sample_name = filename.split('/')[1]
            if sample_name not in md5_sample_list:
                md5_sample_list.append(sample_name)
        md5_sample_list_sorted = sorted(md5_sample_list)
        sample_list_sorted = sorted(sample_list.split(','))
        if md5_sample_list_sorted != sample_list_sorted :
            comm.color_print("交付样本信息与md5文件中的信息不一致")
            sys.exit(1)
        else:
            return True

        
class Path_Deal():
    def __init__( self, release_dir ):
        if not os.path.exists( release_dir ):
            comm.color_print("交付路径不存在，请重新交付")
            sys.exit(1)
        self.release_dir = release_dir
        self.md5 = self.release_dir+'/md5.txt'
        self.raw = self.release_dir+'/Rawdata'
        self.clean = self.release_dir+'/Cleandata'
    
    def check_md5_exists( self ):
        if os.path.exists(self.md5):
            return True
        else:
            return False
    
    def check_raw_clean(self):
        if os.path.exists( self.raw ) and os.path.exists(self.clean ):
            for indir in [ self.raw, self.clean ]:
                self.judge_link(indir)
                self.judge_empty_dir(indir )
                self.judge_empty_file(indir )
            clean_sample_list = self.get_sampleList(self.clean)
            raw_sample_list = self.get_sampleList(self.raw)
            if sorted(clean_sample_list) == sorted(raw_sample_list):
                return len(clean_sample_list), len(raw_sample_list)
            else:
                comm.color_print("Cleandata 和 Rawdata 下的样本不一致，请确认")
                sys.exit(1)
                
        elif os.path.exists( self.raw ) and not os.path.exists(self.clean ) :
            self.judge_link(self.raw)
            self.judge_empty_dir(self.raw )
            self.judge_empty_file(self.raw )
            samplelist = self.get_sampleList(self.raw)
            return 0,len(samplelist)
            
        elif os.path.exists( self.clean ) and not os.path.exists( self.raw ):
            self.judge_link(self.clean)
            self.judge_empty_dir(self.clean )
            self.judge_empty_file(self.clean )
            samplelist = self.get_sampleList(self.clean)
            return len(samplelist),0
        else:
            comm.color_print("交付目录下既没有Cleandata也没有Rawdata，请核实")
            return 0,0
            
    
    def judge_link( self, indir ):
        '''
        判断交付目录/Rawdata/*/* 和交付目录/Cleandata/*/* 的软链接文件是否存在
        返回软连接失效的文件
        '''
        no_link_file = []
        data_file_list = glob.glob(indir +'/*/*')
        for infile in data_file_list :
            if not os.path.exists( infile ):
                no_link_file.append(os.path.basename(infile))
        if no_link_file :
            comm.color_print("{0} 有失效链接，请确认,{1}".format(os.path.basename(indir),no_link_file))
        else:
            comm.color_print(content="{0} 没有失效软链接".format(os.path.basename(indir)),backgroud='40',font_color='37',asc_control='0')

    def judge_empty_dir(self, indir ):
        '''
        判断Rawdata/*和Cleandata/*下面样本目录下是否为空目录
        '''
        empty_list = []
        sample_dir = glob.glob(indir+'/*')
        for sample in sample_dir :
            if os.path.isfile(sample) : continue
            sample_file = glob.glob('{0}/*'.format( sample ))
            if len(sample_file ) == 0 :
                empty_list.append( os.path.basename(sample ))
        if empty_list:
            comm.color_print("{0} 下的样本目录有空目录，请确认,{1}".format(os.path.basename(indir),empty_list))
        else:
            comm.color_print(content="{0} 下没有空目录".format(os.path.basename(indir)),backgroud='40',font_color='37',asc_control='0')
        
    def judge_empty_file(self, indir ):
        empty_list = []
        sample_dir = glob.glob(indir+'/*')
        for sample in sample_dir :
            if os.path.isfile(sample) : continue
            sample_file = glob.glob('{0}/*'.format( sample ))
            for fq_file in sample_file:
                if os.path.getsize(fq_file):continue
                empty_list.append( os.path.basename(sample ))
        
        if empty_list:
            comm.color_print("{0} 下的样本目录有大小为0的文件,请确认,{1}".format(os.path.basename(indir), empty_list))
            sys.exit(1)   
        else:
            comm.color_print(content="{0} 下没有大小为0的文件".format(os.path.basename(indir)),backgroud='40',font_color='37',asc_control='0')
    def get_sampleList( self, indir ):
        '''
        获取样本名，如果clean和raw目录下的样本名不一致，报错退出
        如果不存在clean也不存在raw目录，报错退出
        '''
        sampleDir_list = glob.glob(indir+'/*')
        sample_list = []
        if len(sampleDir_list) < 1 :
            comm.color_print("{0} 没有样本信息".format( indir ))
            return ['']
        for sampleDir in sampleDir_list:
            if os.path.isfile(sampleDir) : continue
            sample_name = os.path.basename(sampleDir)
            if sample_name not in sample_list:
                sample_list.append( sample_name )
        return sample_list
def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument( '-d', '--dir', help='dir to check', dest = 'dir',default=os.getcwd())
    args = parser.parse_args()

    if not args.dir :
        check_dir = os.getcwd()
    else:
        check_dir = os.path.abspath(args.dir)
    
    comm.color_print(content="**正在复核的路径为 {0}\n".format( check_dir ), backgroud='40',font_color='37',asc_control='0')
    my_path_job = Path_Deal( check_dir )
    clean_num, raw_num = my_path_job.check_raw_clean()
    comm.color_print(content="Cleandata下样本数量为 {0}\nRawdata下样本数量为{1}".format( clean_num, raw_num ),backgroud='40',font_color='37',asc_control='0') 
    my_path_job.check_md5_exists()


    my_md5_job = MD5(check_dir)
    my_md5_job.check_repeat_md5()


if __name__ == '__main__':
    main()