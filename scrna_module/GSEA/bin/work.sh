#GSEA分析-输出rds和结果表
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9 /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/gsva.sif Rscript /work/share/acuhtwkcu9/tuchengfang/05_tool/public/GSEA/bin/gsea_msigdb.r -i 差异基因文件 -o 输出目录 --prefix 输出前缀 -t SYMBOL -g human/mouse
-i : 输入文件，第一列必须为gene列，必须有Log2FC列
-g :物种，human /mouse
#GSEA绘图-输出经典图
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9 /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/gsva.sif Rscript /work/share/acuhtwkcu9/tuchengfang/05_tool/public/GSEA/bin/gsea_plot.r -r rds_file -i ids_file -o 输出目录 --prefix 输出结果前缀 --header False
-id文件不包含表头，包含需要绘图的ID信息

#gsea_plot_v2.r 
可以将pd（p值等）使用外部文件传入
