library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'rds', 'r', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'prefix' , 'n', 1, 'character', '输出的文件名'
),
  byrow=T,ncol=5
)


## 读取参数
args=getopt(command)

print_usage <- function(para=NULL){
    cat("参数说明：
    -r:rds文件
    -i:cluster和细胞的类型的对应关系表 cluster,cell_type,sub_type,notes
    -o:输出文件

    使用实例：
    /public/software/apps/singularity/3.7.3/bin/singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript plot_cell_stack.r -r rds -i releation.csv -o rename.rds --infile:cluster和细胞的类型的对应关系表 cluster,cell_type,sub_type,notes

    ")
}

if ( !is.null(args$help) )	{ print_usage(para) }

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$prefix)) {
  cat("Usage: Rscript cluster_umap -i input.rds -o outdir -p prefix\n")
  cat(paste(getopt(command, usage = T), "\n"))
    q() 
}

if (!file.exists(args$input)) {
  cat("Error: input file not exists!\n")
  q()
}

if (!dir.exists(args$outdir)) {
  cat("Error: outdir not exists!\n")
  dir.create(args$outdir, showWarnings = FALSE , recursive = TRUE)
}


library(qusage)
library(scCustomize)
library(patchwork)

out_pre <- paste(args$outdir, args$prefix, sep="/")


pbmc_t <- readRDS( args$rds )
annotation_info <- Pull_Cluster_Annotation(annotation = args$infile)
pbmc <-  Rename_Clusters(seurat_object = pbmc_t, new_idents = annotation_info$new_cluster_idents)
pbmc@meta.data$CellType <- Idents(pbmc)

rds_file <- paste(out_pre, ".rename.rds", sep="")
saveRDS( pbmc, rds_file)

#绘图
pdf_file <- paste(out_pre, ".dimplot.pdf", sep="")

a<- DiscretePalette_scCustomize(num_colors = 26, palette = "alphabet")
#a[3] <- "#E4E1E3FF"
p1 <- DimPlot_scCustom(pbmc_new, reduction = "umap", label = TRUE , colors_use =a)
#p2 <- DimPlot_scCustom(pbmc_new, reduction = "tsne", label = TRUE , colors_use = a)

p3 <- DimPlot_scCustom(pbmc_new, reduction = "umap", label = TRUE , group.by = "seurat_clusters" , colors_use = a)
#p4 <- DimPlot_scCustom(pbmc_new, reduction = "tsne", label = TRUE , group.by = "seurat_clusters" , colors_use = a)

#pdf("cluster_dimplot.pdf", width = 20, height = 10)
p <- p1
print(p)
p <- p3
print(p)
dev.off()

