### 模块： mk_kegg

*模块功能：GO和KEGG富集分析
*模块版本：v1.0.0
*邮箱： chengfangtu@genome.cn

### 使用示例及参数说明：

Usage:
	make -f ${file} indir= outdir= gene_list= prefix= GetGoList GO
	功能说明：
		对给定的gene_list做GO富集分析
	参数说明：
		indir: 含有Gene_All.GO.xls的文件的目录
		outdir:输出目录
		gene_list:待做GO的基因列表,表头为Gene
		prefix:结果文件前缀名
	
	make -f ${file} indir= outdir= gene_list= category= prefix= GetKEGGList KEGG
	功能说明：
		对给定的gene_list做KEGG富集分析
	参数说明：
		indir: 含有Gene_All.GO.xls的文件的目录
		outdir:输出目录
		gene_list:待做GO的基因列表,表头为Gene
		category:物种类型，[fungi,plant,animal]
		prefix:结果文件前缀名

### 输入文件示例
见test/input/
.
├── Gene_All.GO.xls		   基因和GO的对应关系
├── Gene_All.KEGG.xls	   基因和KEGG的对应关系
|-- gene.list              待分析的基因列表

### 运行环境及软件：
	北京238 R python3

### 资源消耗及运行时长
	申请CPU：1
	申请内存：1G
	运行时长：5min

### 输出文件示例
会在输出目录生成GO和KEGG目录，可以直接交付GO和KEGG目录
.
├── GO
│   ├── readme.doc
│   ├── test.go.barplot.pdf
│   ├── test.go.dotplot.pdf
│   ├── test.go.enrichment.xls
│   └── test.go.report.xls
├── go.list
├── KEGG
│   ├── readme.doc
│   ├── test.kegg.barplot.pdf
│   ├── test.kegg.dotplot.pdf
│   ├── test.kegg.enrichment.xls
│   └── test.kegg.report.xls
└── ko.list

GO主要结果文件说明：
1. GO/*.go.report.xls : 候选基因GO统计结果
（1）ID：GO Term的ID
（2）Ontology：该Term 所属分类
（3）Description：GO Term的描述
（4）Count1：富集到该Term的基因数目；
（5）Count2：用于富集分析的基因数目；
（6）Count3：富集到该Term的的背景基因数目；
（7）Count4：用于富集分析时的背景基因数目；
（8）pval：检验后的p值；
（9）p.adjust：BH方法校正后的p值；
（10）qval：检验后的q值；
（11）*Gene：富集到该Term上的基因
（12）*Count：富集到该Term上的基因数目
（13）Links：该GO Term的数据库链接；
（14）Result：该Term是否显著富集，yes，为显著；no，为不显著。

2. GO/*.dotplot.p* : GO富集分析气泡图
选取每个类别最显著的10个GO条目（如果不足10则用该类别全部条目）用气泡图展示
纵坐标表示GO条目，横坐标表示富集到该条目的基因数量占总基因的比例，颜色表示padjust，颜色越红表示越显著；气泡大小表示富集到该条目的基因数量，气泡越大表示基因数量越多。

3. GO/*.barplot.p* : GO富集分析气泡
选取每个类别最显著的10个GO条目（如果不足10则用该类别全部条目）用条形图展示
纵坐标表示GO条目，横坐标表示富集到该条目的基因数量，颜色表示padjust，颜色越红表示越显著。

KEGG主要结果文件说明：
1. KEGG/*.kegg.report.xls : 候选基因KEGG统计结果表
（1）Map：kegg通路编号
（2）Name：kegg通路名称
（3）Count1：富集到该通路的基因数目；
（4）Count2：用于富集分析的基因数目；
（5）Count3：富集到该通路的的背景基因数目；
（6）Count4：用于富集分析时的背景基因数目；
（7）pval：检验后的p值；
（8）p.adjust：BH方法校正后的p值；
（8）qval：检验后的q值；
（9）*Gene：富集到该通路上的基因
（10）*Count：富集到该通路上的基因数目
（11）Links：该map的数据库链接；
（12）Result：该map是否显著富集，yes，为显著；no，为不显著。

2. KEGG/*.dotplot.p* : KEGG富集分析气泡图
选取最显著的30个pathway通路（如果不足30则用全部通路）用气泡图展示
纵坐标表示通路名称，横坐标表示富集到该通路的基因数量占总基因的比例，颜色表示padjust，颜色越红表示越显著；气泡大小表示富集到该通路的基因数量，气泡越大表示基因数量越多。

3. KEGG/*.barplot.p* : KEGG富集分析条形图
选取最显著的30个pathway通路（如果不足30则用全部通路）用条形图展示
纵坐标表示通路名称，横坐标表示富集到该通路的基因数量，颜色表示padjust，颜色越红表示越显著。




### 注意事项
投递的时候有问题。