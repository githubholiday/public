file=$(abspath $(firstword $(MAKEFILE_LIST)))
BIN=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))/
ifeq ($(strip $(config)),)
	Bconfig=$(BIN)/../config/config.txt
else
	Bconfig=$(config)
endif
include $(Bconfig)

Help:
	@echo Description:
	@echo -e "\t" 该脚本用于做GO和KEGG的富集分析
	@echo Usage:
	@echo -e "\t" make -f ${file} indir= outdir= gene_list= prefix= GetGoList GO
	@echo Parameters:
	@echo -e "\t" indir: 含有Gene_All.GO.xls的文件的目录
	@echo -e "\t" outdir:输出目录
	@echo -e "\t" gene_list:待做GO的基因列表,表头为Gene
	@echo -e "\t" prefix:结果文件前缀名
	
	@echo -e "\n" Usage:
	@echo -e "\t" make -f ${file} indir= outdir= gene_list= category= prefix= GetKEGGList KEGG
	@echo Parameters:
	@echo -e "\t" indir: 含有Gene_All.GO.xls的文件的目录
	@echo -e "\t" outdir:输出目录
	@echo -e "\t" gene_list:待做GO的基因列表,表头为Gene
	@echo -e "\t" category:物种类型，[fungi,plant,animal]
	@echo -e "\t" prefix:结果文件前缀名

infile=$(indir)/Gene_All.GO.xls
.PHONY:GetGoList
GetGoList:
	echo "########## GetGoList start at" `date`
	mkdir -p $(outdir)
	sed '1d' $(infile) > $(outdir)/go.list
	echo "########## GetGoList end at" `date`

go_dir=$(outdir)/GO
go=$(outdir)/go.list
.PHONY:GO
GO:
	echo "########## GO Clusterprofiler start at" `date`
	mkdir -p $(go_dir)
	$(SinRun) $(function_sif) make -f $(BIN)/GO_clusterProfiler/GO_clusterProfiler.mk species=no term2gene=$(go) prefix=$(prefix) outdir=$(go_dir) genelist=$(gene_list) config=$(Bconfig) GO_clusterProfiler
	rm -r $(go_dir)/*result
	rm -r $(go_dir)/*example*
	rm -r $(go_dir)/*list
	cp -r $(BIN)/../doc/go.readme.doc $(go_dir)/readme.doc
	echo "########## GO Clusterprofiler end at" `date`

infile=$(indir)/Gene_All.KEGG.xls
.PHONY:GetKEGGList
GetKEGGList:
	echo "########## GetKEGGList start at" `date`
	mkdir -p $(outdir)
	sed '1d' $(infile) > $(outdir)/ko.list
	echo "########## GetKEGGList end at" `date`

kegg_dir=$(outdir)/KEGG
ko=$(outdir)/ko.list
.PHONY:KEGG
KEGG:
	echo "########## KEGG Clusterprofiler start at" `date`
	mkdir -p $(kegg_dir)
	$(SinRun) $(function_sif) make -f $(BIN)/KEGG_clusterProfiler/KEGG_clusterProfiler.mk species=no term2gene=$(ko) prefix=$(prefix) outdir=$(kegg_dir) genelist=$(gene_list) config=$(Bconfig) category=$(category) KEGG_clusterProfiler
	rm -r $(kegg_dir)/*result
	rm -r $(kegg_dir)/*id
	rm -r $(kegg_dir)/*example*
	rm -r $(kegg_dir)/*list
	cp -r $(BIN)/../doc/kegg.readme.doc $(kegg_dir)/readme.doc
	echo "########## KEGG Clusterprofiler end at" `date`