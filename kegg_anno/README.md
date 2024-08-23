### 模块： mk_RDA_CCA

*模块功能：KEGG注释结果整理、统计、汇总
*模块版本：v1.0.0
*邮箱： chengfangtu@genome.cn

### 使用示例及参数说明：

Usage:
    make -f $(makefile_name) input= kegg_level_file= TPM= outdir= KEGG_Combine
	功能说明：
	    对切割的fa做的kegg注释的结果进行合并以及增加TPM值
	参数说明：
	    config: [文件|可选]  模块配置文件，和软件相关参数，默认为$(makefile_dir)/config/config.txt 
	    input: [文件|必需]  KEGG注释结果,一般为*.out
	    kegg_level_file: [路径|必需]  KEGG 不同Level对应关系表，一般为/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/Database/KEGG/current/data/map_pathway.list
	    TPM: [字符|必需]  基因在不同样本中的TPM表，第一行为表头，第一列为gene名，值为TPM值
	    outdir: [字符|必需] 结果输出目录

    make -f $(makefile_name) infile= outdir= KEGG_level_tpm_stat
	功能说明：
	    对KEGG注释后的不同Level计算丰度以及输出基因列表
	参数说明：
	    config: [文件|可选]  模块配置文件，和软件相关参数，默认为$(makefile_dir)/config/config.txt 
	    infile: [文件|必需]  KEGG注释结果，也就是KEGG_Combine的输出结果
	    outdir: [路径|必需]  结果输出路径，会输出不同level的丰度和对应的基因列表

### 输入文件示例
见test/input/
.
├── *.out           KEGG注释后的文件
├── tpm.xls         所有样本的TPM表格

### 运行环境及软件：
	北京238 R python3

### 资源消耗及运行时长
	申请CPU：1
	申请内存：1G
	运行时长：5min

### 输出文件示例
.
├── RDA_CCA.coordinate.env.xls       RDA/CCA结果图中环境因子的坐标 
├── RDA_CCA.coordinate.pdf/png       RDA/CCA结果图，根据DCA的结果进行选择，只出一种图，会在图里标明是RDA或者CCA
└── RDA_CCA.coordinate.sample.xls    RDA/CCA结果图中样本的坐标

主要结果文件说明：
（1）RDA_CCA.coordinate.pdf/png
1）环境向量的长度表示样方物种的分布与该环境因子相关性的大小，长度越长，相关性越大；
2）环境向量与约束轴夹角的大小表示环境因子与约束轴相关性的大小，夹角小说明关系密切，若正交则不相关；
3）样本点与箭头距离越近，该环境因子对样本的作用越强；
4）样本位于箭头同方向，表示环境因子与样本物种群落的变化正相关，样本位于箭头的反方向，表示环境因子与样本物种群落的变化负相关。

（2）RDA_CCA.coordinate.env.xls
CCA1/RDA1: 第一轴的坐标
CCA2/RDA2: 第二轴的坐标
factor: 环境因子

（3）RDA_CCA.coordinate.sample.xls 
CCA1/RDA1: 第一轴的坐标
CCA2/RDA2: 第二轴的坐标
Sample: 样本 
Group: 分组

### 注意事项
投递的时候有问题。