### 模块： mk_kegg

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
├── TPM.example.xls         所有样本的TPM表格

### 运行环境及软件：
	北京238 R python3

### 资源消耗及运行时长
	申请CPU：1
	申请内存：1G
	运行时长：5min

### 输出文件示例
会在输出目录生成kegg_upload目录，将最终要交付的文件进行整理
.
├── Gene_All.KEGG.Level1.genelist.xls 
├── Gene_All.KEGG.Level1.TPM.Summary.xls 
├── Gene_All.KEGG.Level2.genelist.xls 
├── Gene_All.KEGG.Level2.TPM.Summary.xls 
├── Gene_All.KEGG.Level3.genelist.xls 
├── Gene_All.KEGG.Level3.TPM.Summary.xls 
├── Gene_All.KEGG.TPM.xls 
└── readme.doc

主要结果文件说明：
（1）Gene_All.KEGG.Level*.genelist.xls
KEGG不同层级的上注释到的基因信息，第一列为KEGG层级，最后一列为gene_id（基因之间以|分割）
如果是Level3则会将其对应的Level2和Level1也标注在表格中
 (2)Gene_All.KEGG.Level1.TPM.Summary.xls 
 KEGG不同层级在不同样本中的TPM值（该条目上所有基因的TPM值总和）

（3）Gene_All.KEGG.TPM.xls
基因注释到KEGG数据库的结果总表
Gene_ID:基因ID序号
Annotation:基因注释到的K编号
Map:基因注释到的Map编号
Level3:基因注释到KEGG数据库的Level3层级（即pathway）
Level2:基因注释到KEGG数据库的Level2层级
Level1:基因注释到KEGG数据库的Level1层级
样本：后面所有列为该样本中基因的TPM值

### 注意事项
投递的时候有问题。