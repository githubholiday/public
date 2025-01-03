#species?=hg19
tmpdir=$(dir  $(abspath $(firstword $(MAKEFILE_LIST))))
BIN=$(tmpdir)../
bin=$(tmpdir)
include $(tmpdir)/../config/config.txt
include $(ref_index)/$(species).txt

comma:=,
name=$(subst $(comma), ,$(compare))
cmp=$(subst $(comma),_,$(compare))
de_report=$(indir)/cmp/$(tool)/$(cmp)//$(cmp).report.xls
kegg_dir=$(outdir)/Function/KEGG/$(cmp)
de = $(kegg_dir)/de.list
# DrawKEGGHeat
heat_dir=$(outdir)/Function/KEGG/heat
kegg_files=$(outdir)/Function/KEGG/*/*kegg.report.xls
select=15
#KEGGMap
#de_report
output=$(kegg_dir)/$(cmp)
up=$(kegg_dir)/up.list
down=$(kegg_dir)/down.list
mapdb=$(KEGG_Base)/kegg_map
ko2map=$(KEGG_Base)/ko2map/ko2map.xls

HELP:
	@echo 功能:有参转录组流程进行KEGG富集分析
	@echo 模块说明:
	@echo -e "\t"KEGGCandidate:使用de_repoirt文件准备KEGG_Enrich,KEGGMap的输入文件，如果不使用该模块，后面各个模块输入要求的输入文件即可
	@echo -e "\t"KEGG_Enrich:基于clusterprofiler模块进行KEGG分析，输入为差异分析的report.xls文件
	@echo -e "\t"DrawKEGGHeat:基于多比较组的 *.kegg.report.xls文件绘制多比较组的富集分析热图，默认选取前15行
	@echo -e "\t" DrawKEGGHeat 用法 make -f makefile config= kegg_file=*kegg.report.xls  heatdir= select= DrawKEGGHeat
	@echo -e "\t" KEGGMap 用法1 make -f makefile config= up= down=  kegg_dir= kegg_report= cmp=T_C KEGGMap
	@echo -e "\t" KEGGMap 用法2,通过outdir匹配输出路径,默认kegg_dir=outdir/Function/KEGG/cmp，kegg_report=kegg_dir/cmp.kegg.report.xls
	@echo -e "\t" KEGGECandidate make -f makefile config= de_report= outdir= compare= KEGGCandidate
	@echo -e "\t" KEGGEnrich 用法1 make -f makefile config= de= outdir= compare= KEGGEnrich
	@echo -e "\t" KEGGEnrich 用法2 make -f makefile config= indir= outdir= compare=T,C  tool=deseq2/degseq KEGGCandidate KEGGEnrich
	@echo -e 参数说明：
	@echo -e "\t"config:配置文件，包括软件、数据库、路径等配置，默认为流程目录下的../../software/software.txt
	@echo -e "\t"tool:deseq2或者degseq，主要用于识别de_report文件路径
	@echo -e "\t"config:配置文件，包括软件、数据库、路径等配置，默认为流程目录下的./config.txt
	@echo -e "\t"kegg_file：各比较组的kegg cluster profiler的report.xls文件列表
	@echo -e "\t"heat_dir：比较组热图输出路径
	@echo -e "\t"select：[select]绘制热图是选取的数目，默认是15
	@echo -e "\t"up:[must]显著差异上调基因列表，无表头
	@echo -e "\t"down:[must]显著差异下调基因列表，无表头
	@echo -e "\t"cmp:[must]比较组名，如T_C，也是KEGGMap模块的输出文件的前缀名
	@echo -e "\t"kegg_report:[must]KEGGMap模块中的输入文件，为KEGG clusterprofiler的输出report文件
	@echo -e "\t"kegg_dir:[must]KEGGMap模块中的输出路径
	@echo -e "\t"outdir:[select]输出路径,会在outdir目录下创建Function/KEGG/cmp目录
	
	@echo 'KEGG Analysis'

KEGG_Candidate:
	echo kegg prepare start at `date`
	[ -d $(kegg_dir) ] || mkdir -p $(kegg_dir)
	grep -v -w 'no' $(de_report)  > $(de)
	sed -i 's/Gene_ID/Gene/g' $(de)
	sed -i 's#\tUp\t#\tup\t#g' $(de)
	sed -i 's#\tDown\t#\tdown\t#g' $(de)
	grep -w 'yes' $(de_report) | grep -w 'Up' | cut -f1 > $(up) 
	grep -w 'yes' $(de_report) | grep -w 'Down' | cut -f1 > $(down)
	echo kegg prepare finish at `date`

KEGG_Enrich:
	echo kegg enrich start at `date`
	[ -d $(kegg_dir) ] || mkdir -p $(kegg_dir)
	if [ `wc -l $(de) | awk '{print $$1}'` -gt 1 ];\
		then \
			make -f $(BIN)/KEGG_clusterProfiler/KEGG_clusterProfiler.mk species=no term2gene=$(KEGG_annotate) prefix=$(cmp) outdir=$(kegg_dir) genelist=$(de) category=$(category) config=$(BIN)/config/config.txt KEGG_clusterProfiler ;\
		else \
			echo "no de gene" ;\
	fi;
	echo kegg enrich finish at `date`

DrawKEGGHeat:
	echo KEGG HeatMap start at `date`
	mkdir -p $(heat_dir)
	$(PYTHON3) $(BIN)/KEGG/overlap_kegg.py -f $$(ls $(kegg_files)|while read i;do if  [[ `wc -l $$i|cut -d " " -f 1` -gt 1 ]]; then echo $$i;fi ;done) -i 1 -c 7 -head 1 -t kegg  > $(heat_dir)/keggheat.xls
	if [ `wc -l $(heat_dir)/keggheat.xls | awk '{print $$1}'` -gt 2 ] ;\
	then \
		$(RSCRIPT2) $(BIN)/KEGG/r/pheatmap.r $(heat_dir)/keggheat.xls $(heat_dir)/keggheat.pdf KEGG ;\
		convert $(heat_dir)/keggheat.pdf $(heat_dir)/keggheat.png ;\
		head -$(select) $(heat_dir)/keggheat.xls > $(heat_dir)/keggheat.tmp.xls ;\
		$(RSCRIPT2) $(BIN)/KEGG/r/pheatmap.r $(heat_dir)/keggheat.tmp.xls $(heat_dir)/example.keggheat.pdf KEGG ;\
		convert $(heat_dir)/example.keggheat.pdf $(heat_dir)/example.keggheat.png ;\
	else \
		echo 'none compare group got Significant KEGG map Enrichment!' ;\
		echo [Warning:] No example.keggheat.pdf generated at all! >> $(log_file);\
	fi
	echo KEGG HeatMap end at `date`

KEGGMap:
	[ -s $(dir $(kegg_dir))/kegg_picture ] && rm -r $(dir $(kegg_dir))/kegg_picture || echo dir is clean
	[ -s $(dir $(kegg_dir))/mapcolor ] && rm -r $(dir $(kegg_dir))/mapcolor || echo dir is clean
	mkdir -p $(dir $(kegg_dir))/kegg_picture
	mkdir -p $(dir $(kegg_dir))/mapcolor
	if [ $(denovo) ] ;\
		then \
		$(PERL) $(bin)/keggFileConfig.pl -i $(anno_file) -o $(kegg_dir)/$(cmp).kegg_file_config ;\
	else \
		$(PERL) $(bin)/keggFileConfig_new.pl -up $(up) -down $(down) -ko2map $(ko2map) -keg $(KEGG_Base)/ko2map/ko00001.keg -ko $(KEGG_annotate) -out $(prefix).kegg_file_config ;\
	fi
	$(PYTHON3) $(bin)/HtmlTable.py -i $(kegg_report) -o $(prefix).html_tmp -c 2
	$(PERL) $(bin)/map_color.pl -kegg_report $(kegg_report) -mapdb $(mapdb) -od $(dir $(kegg_dir))/mapcolor
	$(PYTHON3) $(bin)/kegg.py -r $(kegg_report) -f $(kegg_dir).kegg_file_config -m $(dir $(kegg_dir))/mapcolor -o $(dir $(kegg_dir))/kegg_picture
	cp -r $(bin)/kegg_html/src $(dir $(kegg_dir))/

