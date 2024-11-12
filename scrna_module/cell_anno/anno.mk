
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

subset_dir=$(outdir)/subset
Idents=CellType
.PHONY:merge
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
marker_dir=cluster_dir=$(outdir)/marker
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
	$(SinRun) $(SIF) Rscript $(ScriptDir)/rename_celltype.r -r $(inrds) -i $(rename_file) -o $(rename_dir) -p $(prefix) -O oldidents -N newidents
	echo "########## rename rds end at" `date`

plot_dir=$(outdir)/draw
maker_plot:
	$(SinRun) $(SIF) Rscript $(ScriptDir)/draw_plot.r -i $(inrds) -o $(plot_dir) -n $(prefix)  -g $(gmt)

## B细分
### 提取对应细胞
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script/subset_rds.r -i ../new.rds -o ./ -n B -c B -I CellType

-c:
### 对提取后rds cluster
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script/subset_rds_cluster.r -i B.subset.rds  -o . -n B
### findallmarker 
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script/findallmarker.r -i B_0.6.umap.rds -o . -n B -I seurat_clusters
### 利用marker画图 

singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script//draw_plot.r -i B_0.6.umap.rds -o . -n B  -g gene.gmt 

### 画类似于python的cluster dotplot图 
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/sceasy.sif Rscript /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script/toAdata.r -i B_0.6.umap.rds -o . -n B_0.6.umap.adata
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/scanpy/scanpy.sif python  /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script/draw_cluster_dotplot.py -i B_0.6.umap.adata.h5ad -g  gene.gmt -o B_all_gene_dotplot_cluster.pdf
 

###rename 
### 修改配置文件，填入细胞类型 
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript /work/share/acuhtwkcu9/taoxiao/05_tool/my_module/public/scrna_module/cell_anno/script/rename_celltype.r -r rds -i cell.csv -o outdir -p prefix -O oldidents -N newidents

### replace name back to original rds

singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script/replace_name_from_small_to_big.r -s B_0.6.umap.rds -o . -n B -c  seurat_clusters  -b  ../new.rds 

# 如果是subset，
# 提取subset,cluster
# 1 findallmarker+function
# 2 marker绘图 =>确认分群
# 3 rename+绘图
# 4 findallmarker
# 5 function

首次：
findallmarker+function+marker绘图

subset:
提取subset+cluster+findallmarker+function+maker绘图+