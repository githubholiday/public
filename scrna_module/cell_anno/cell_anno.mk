
file=$(abspath $(firstword $(MAKEFILE_LIST)))
BIN=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))/
ifeq ($(strip $(config)),)
	Bconfig=$(BIN)/config/config.txt
else
	Bconfig=$(config)
endif
include $(Bconfig)
ScriptDir=$(BIN)/script

Help:
	@echo Description:
	@echo -e "\t" 单细胞转录组的基础分析脚本
	@echo Usage:
	@echo -e "\n" Usage:合并rds或者三个matrix文件
	@echo -e "\t" make -f ${file} sample_info= outdir= sample= merge
	@echo Parameters:
	@echo -e "\t" sample_info: 样本信息表，第一行为name\tpath\tgroup
	@echo -e "\t" sample: 样本名
	@echo -e "\t" outdir: 输出目录
	@echo -e "\t" config: 配置文件

	@echo -e "\n" Usage:对合并后的数据进行过滤、
	@echo -e "\t" make -f ${file} inrds= sample= outdir= qc
	@echo Parameters:
	@echo -e "\t" inrds: 合并样本的rds文件
	@echo -e "\t" sample: 样本名称
	@echo -e "\t" outdir: 输出目录，并在输出目录下创建QC 

	@echo -e "\n" Usage:对合并后的数据进行聚类
	@echo -e "\t" make -f ${file} inrds= sample= outdir= resolution= cluster
	@echo Parameters:
	@echo -e "\t" inrds:rds文件，一般为qc后的rds
	@echo -e "\t" sample:样本名称
	@echo -e "\t" outdir:输出目录，并在输出目录下创建 cluster_resolution
	@echo -e "\t" resolution:分辨率

subset_dir=$(outdir)/subset
Idents=CellType
.PHONY:subset
subset:
	echo "########## subset start at" `date`
	mkdir $(merge_dir)
	$(SinRun) $(SIF) Rscript $(ScriptDir)/subset_rds -i $(inrds) -o $(subset_dir) -n $(prefix) -c $(cluster) -I $(Idents)
	echo "########## subset end at" `date`

inrds=$(outdir)/subset/$(prefix).subset.rds
cluster_dir=$(outdir)/cluster
subset_cluster:
	echo "########## subset cluster start at" `date`
	mkdir -p $(cluster_dir)
	$(SinRun) $(SIF) Rscript $(ScriptDir)/subset_rds_cluster.r -i $(inrds) -o $(cluster_dir) -n $(prefix)
	echo "########## subset cluster  end at" `date`

inrds=$(outdir)/subset/$(prefix).subset.rds
marker_dir=$(outdir)/de
Idents=seurat_clusters
findallmarker:
	echo "########## findallmarker start at" `date`
	mkdir -p $(marker_dir)
	$(SinRun) $(SIF) Rscript $(ScriptDir)/findallmarker.r -i $(inrds) -o $(marker_dir) -n $(prefix) -I $(Idents)
	echo "########## findallmarker end at" `date`

rename_dir=$(outdir)/rename
oldidents=seurat_clusters
newidents=CellType
rename:
	echo "########## rename rds start at" `date`
	mkdir -p $(rename_dir)
	$(SinRun) $(SIF) Rscript $(ScriptDir)/rename_celltype.r -r $(inrds) -i $(rename_file) -o $(rename_dir) -p $(prefix) -O $(oldidents) -N $(newidents)
	echo "########## rename rds end at" `date`

plot_dir=$(outdir)/marker_draw
marker_plot:
	echo "########## marker plot start at" `date`
	mkdir -p $(plot_dir)
	$(SinRun) $(SIF) Rscript $(ScriptDir)/draw_plot.r -i $(inrds) -o $(plot_dir) -n $(prefix)  -g $(gmt) -a $(assays)
	echo "########## marker plot end at" `date`


plot_dir=$(outdir)/marker_draw
maker_plot_py:
	echo "########## marker plot start at" `date`
	mkdir -p $(plot_dir)
	$(SinRun) $(SIF2) Rscript $(ScriptDir)/toAdata.r -i $(inrds) -o $(plot_dir) -n $(prefix).adata
	$(SinRun) $(scanpy_sif) python  $(ScriptDir)/draw_cluster_dotplot.py -i $(plot_dir)/$(prefix).adata.h5ad -g $(gmt) -o $(outdir)/$(prefix).all_gene.dotplot.pdf
	echo "########## marker plot end at" `date`

### replace name back to original rds
replace:
	singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script/replace_name_from_small_to_big.r -s B_0.6.umap.rds -o . -n B -c  seurat_clusters  -b  ../new.rds 
all=yes
cmp_de_dir=$(outdir)/cmp_de/
.PHONY:cmp_de
cmp_de:
	echo "########## cmp de start at" `date`
	mkdir -p $(cmp_de_dir)/$(prefix)
	mkdir -p $(outdir)/shell/$(prefix)
	$(SinRun) $(SIF) Rscript $(ScriptDir)/deseq_bulk_and_scRNA.r -i $(inrds) -o $(cmp_de_dir)/$(prefix) -n $(prefix) -c $(cmp) -m scRNA -I $(Idents)
	make -f $(BIN)/../Function/v2.0.0/bin/anno.mk infile=$(cmp_de_dir)/$(prefix)/*.xls outdir=$(cmp_de_dir)/$(prefix)/function shell_dir=$(outdir)/shell/$(prefix) species_conf=$(species_conf) all=$(all) Serial_Function
	echo "########## cmp de end at" `date`
all=yes
de_dir=$(outdir)/cluster_de
.PHONY:cluster_de
cluster_de:
	echo "########## cmp de start at" `date`
	mkdir -p $(outdir)/cluster_de
	mkdir -p $(outdir)/shell/
	$(SinRun) $(SIF) Rscript $(ScriptDir)/findallmarker.r -i $(inrds) -o $(de_dir) -n $(prefix) -I $(Idents)
	make -f $(BIN)/../Function/v2.0.0/bin/anno.mk infile=$(de_dir)/*.xls outdir=$(de_dir)function shell_dir=$(outdir)/shell/ species_conf=$(species_conf) all=$(all) Serial_Function
	echo "########## cmp de end at" `date`