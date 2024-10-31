library(getopt)

command=matrix(c( 
    'help', 'h', 0,'logical', '帮助文档',
    'rds', 'r', 1, 'character', '输入的rds文件',
    'input', 'i', 1, 'character', '细胞类型对应关系',
    'outdir', 'o', 1, 'character', '输出的目录',
    'prefix' , 'p', 1, 'character', '输出的文件名',
    'OldIdents' , 'O', 1, 'character', '重命名前的slots名,默认为seurat_clusters',
    'NewIdents' , 'N', 1, 'character', '重命名后的slots名,默认为CellType'
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
    -p:输出的前缀
    -O:重命名前的slot名称
    -I:重命名后的slot名称

    使用实例：
    /public/software/apps/singularity/3.7.3/bin/singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript plot_cell_stack.r -r rds -i releation.csv -o outdir -p prefix -O seurat_clusters -N CellType
    ")
}

if ( !is.null(args$help) )	{ print_usage(para) }

if (!is.null(args$help)|| is.null(args$rds) || is.null(args$input) || is.null(args$outdir) || is.null(args$prefix)) {
  cat("Usage: Rscript cluster_umap -r input.rds -o outdir -p prefix -i input -OI seurat_clusters -NI CellType\n")
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

if (!is.null(args$OldIdents)) {
    old_idents <- args$Idents 
}else{
    old_idents <- "seurat_clusters"
}

if (!is.null(args$NewIdents)) {
    new_idents <- args$NewIdents 
}else{
    new_idents <- "CellType"
}



library(qusage)
library(scCustomize)
library(patchwork)

out_pre <- paste(args$outdir, args$prefix, sep="/")


pbmc_t <- readRDS( args$rds )
#Idents(pbmc_t) <- old_idents

annotation_info <- Pull_Cluster_Annotation(annotation = args$infile)
pbmc <-  Rename_Clusters(seurat_object = pbmc_t, new_idents = annotation_info$new_cluster_idents)

#如果更换名称的话，可以改CellType
pbmc[[new_idents]] <- Idents(pbmc)

rds_file <- paste(out_pre, ".rename.rds", sep="")
saveRDS( pbmc, rds_file)

#绘图
pdf_file <- paste(out_pre, ".dimplot.pdf", sep="")

a<- DiscretePalette_scCustomize(num_colors = 26, palette = "alphabet")
#a[3] <- "#E4E1E3FF"
p1 <- DimPlot_scCustom(pbmc, reduction = "umap", label = TRUE , colors_use =a)
#p2 <- DimPlot_scCustom(pbmc_new, reduction = "tsne", label = TRUE , colors_use = a)

p2 <- DimPlot_scCustom(pbmc, reduction = "umap", label = TRUE , group.by = "seurat_clusters" , colors_use = a)
#p4 <- DimPlot_scCustom(pbmc_new, reduction = "tsne", label = TRUE , group.by = "seurat_clusters" , colors_use = a)
print(p1)
print(p2)
dev.off()

