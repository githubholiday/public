a## B细分
### 提取对应细胞
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script/subset_rds.r -i ../new.rds -o ./ -n B -c B -I CellType

-c:
### 对提取后rds cluster
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script/subset_rds_cluster.r -i B.subset.rds  -o . -n B
### findallmarker 
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script/findallmarker.r -i B_0.6.umap.rds -o . -n B -I seurat_clusters
### 利用marker画图 

singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/script//draw_plot.r -i B_0.6.umap.rds -o . -n B  -g gene.gmt -s F

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