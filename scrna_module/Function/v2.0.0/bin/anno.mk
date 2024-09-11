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
	@echo -e "\t" 该脚本用于三代数据整理交付目录
	@echo Usage:
	@echo -e "\t" make -f ${file} indir= sample_list= release_dir= Releasd_dir
	 @echo  Parameters:
	@echo -e "\t" indir: 过滤目录Filter_Result
	@echo -e "\t" sample_list：样本列表文件
	@echo -e "\t" release_dir：整理输出目录

outfile=$(outdir)/$(sample).xls
annoout=$(outdir)/$(sample).anno.xls
.PHONY:Anno
Anno:
	echo "########## Anno start at" `date`
	mkdir -p $(outdir)
	python3 $(BIN)/anno/add_geneId.py -i $(infile) -r $(relation_file) -o $(outfile)
	$(SinRun) $(anno_sif) perl $(BIN)/anno/add_GO_KEGG.pl -in $(outfile)  -out  $(annoout) -anno $(anno_file) -title NR:Seq-id NR:Score NR:Evalue NR:Description NT:Seq-id NT:Score NT:Evalue NT:Description Uniprot:UniProtKB-AC Uniprot:Score Uniprot:Evalue Uniprot:Description COG:gene COG:Score COG:Eval COG:num Pfam:pfam_ID Pfam:pfam_Name Pfam:pfam_Description GO:biological_process GO:cellular_component GO:molecular_function KEGG:KO KEGG:Description
	echo "########## Anno end at" `date`
#获取上/下调基因
gene_list=$(outdir)/$(up_or_down).gene.xls
all_list=$(outdir)/all.gene.xls
.PHONY:GetList
GetList:
	echo "########## Anno start at" `date`
	mkdir -p $(outdir)
	grep "${up_or_down}"  $(infile) | cut -f 1 > $(gene_list)
	/usr/bin/sed -i '1i Gene' $(gene_list)
	echo "########## Anno end at" `date`

go_dir=$(outdir)/GO
GO:
	echo "########## GO Clusterprofiler start at" `date`
	mkdir -p $(go_dir)
	$(SinRun) $(function_sif) make -f $(BIN)/GO_clusterProfiler/GO_clusterProfiler.mk species=no term2gene=$(go) prefix=$(sample) outdir=$(go_dir) genelist=$(gene_list) config=$(Bconfig) GO_clusterProfiler
	$(PYTHON3) $(BIN)/anno/add_geneName.py -i $(go_dir)/$(sample).go.report.xls -r $(relation_file) -cn Gene -o $(go_dir)/$(sample).go.report.name.xls
	$(PYTHON3) $(BIN)/anno/add_geneName.py -i $(go_dir)/$(sample).go.enrichment.xls -r $(relation_file) -cn geneID -o $(go_dir)/$(sample).go.enrichment.name.xls
	mv $(go_dir)/$(sample).go.enrichment.name.xls $(go_dir)/$(sample).go.enrichment.xls
	mv $(go_dir)/$(sample).go.report.name.xls $(go_dir)/$(sample).go.report.xls
	rm -r $(go_dir)/*result
	rm -r $(go_dir)/*example*
	rm -r $(go_dir)/*list
	cp -r $(BIN)/../report/readme.docx $(go_dir)/../
	echo "########## GO Clusterprofiler end at" `date`

kegg_dir=$(outdir)/KEGG
KEGG:
	echo "########## KEGG Clusterprofiler start at" `date`
	mkdir -p $(kegg_dir)
	$(SinRun) $(function_sif) make -f $(BIN)/KEGG_clusterProfiler/KEGG_clusterProfiler.mk species=no term2gene=$(ko) prefix=$(sample) outdir=$(kegg_dir) genelist=$(gene_list) config=$(Bconfig) category=$(category) KEGG_clusterProfiler
	$(PYTHON3) $(BIN)/anno/add_geneName.py -i $(kegg_dir)/$(sample).kegg.report.xls -r $(relation_file) -cn Gene -o $(kegg_dir)/$(sample).kegg.report.name.xls
	$(PYTHON3) $(BIN)/anno/add_geneName.py -i $(kegg_dir)/$(sample).kegg.enrichment.xls -r $(relation_file) -cn geneID -o $(kegg_dir)/$(sample).kegg.enrichment.name.xls
	mv $(kegg_dir)/$(sample).kegg.enrichment.name.xls $(kegg_dir)/$(sample).kegg.enrichment.xls
	mv $(kegg_dir)/$(sample).kegg.report.name.xls $(kegg_dir)/$(sample).kegg.report.xls
	rm -r $(kegg_dir)/*result
	rm -r $(kegg_dir)/*id
	rm -r $(kegg_dir)/*example*
	rm -r $(kegg_dir)/*list
	echo "########## KEGG Clusterprofiler end at" `date`


