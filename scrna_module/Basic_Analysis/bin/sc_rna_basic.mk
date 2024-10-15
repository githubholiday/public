
file=$(abspath $(firstword $(MAKEFILE_LIST)))
BIN=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))/../
ifeq ($(strip $(config)),)
	Bconfig=$(BIN)/config/config.txt
else
	Bconfig=$(config)
endif
include $(Bconfig)
ScriptDir=$(BIN)/bin/script

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

merge_dir=$(outdir)/merge
.PHONY:merge
merge:
	echo "########## merge rds or matrix start at" `date`
	mkdir $(merge_dir)
	$(SinRun) $(SIF) Rscript  $(ScriptDir)/merge_sample.r -i $(sample_info)  -o $(merge_dir) -n $(sample)
	echo "########## merge rds or matrix end at" `date`

inrds=$(outdir)/merge/$(sample).rds
qc_dir=$(outdir)/QC
.PHONY:qc
qc:
	echo "########## scRNA Basic Analysis start at" `date`
	mkdir -p $(qc_dir)
	$(SinRun) $(SIF) Rscript $(ScriptDir)/qc.r -i $(inrds) -m ^Mt- -o $(qc_dir) -n $(sample)
	echo "########## scRNA Basic Analysis end at" `date`

inrds=$(outdir)/QC/$(sample).after_qc.rds
cluster_dir=$(outdir)/cluster_$(resolution)
cluster_draw_dir=$(cluster_dir)/draw
cluster:
	echo "########## scRNA Basic Analysis start at" `date`
	mkdir -p $(cluster_draw_dir)
	$(SinRun) $(SIF) Rscript $(ScriptDir)/qc.r -i $(inrds) -m ^Mt- -o $(cluster_dir) -n rat_merge
	$(SinRun) $(SIF) Rscript $(ScriptDir)/cluster_umap.r -i $(inrds) -o $(cluster_dir) -n $(sample)  -r $(resolution) 
	#$(SinRun) /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/louper.sif Rscript $(ScriptDir)/toloupr.r -i $(cluster_dir)/$(sample).rds -o $(cluster_dir) -n $(sample).cloupe
	#$(SinRun) $(SIF) Rscript $(ScriptDir)/findallmarker.r -i $(cluster_dir)/$(sample).rds -o $(cluster_dir) -n $(sample)
	#$(SinRun) $(SIF) Rscript $(ScriptDir)/draw_cluster_composition.r -i $(cluster_dir)/$(sample).rds -o $(cluster_draw_dir)  -n $(sample)
	echo "########## scRNA Basic Analysis end at" `date`

marker_plot:

	$(SinRun) $(SIF) Rscript $(ScriptDir)/draw_plot.r -i $(cluster_dir)/$(sample).rds -o $(cluster_draw_dir)  -n $(sample) -g $(gene_marker) -r umap.harmony
	########## scRNA Basic Analysis end at" `date`
