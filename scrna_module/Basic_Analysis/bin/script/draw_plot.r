library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名',
  'gmt' , 'g', 1, 'character', 'gmt文件',
  'assay', 'a', 1, 'character', 'assay 标签名，可选RNA Spatial SCT',
  'reduction' , 'r', 2, 'character', '降维的方法，默认为"umap", 可选为"umap"、"tsne" 、"integrated.cca"、"integrated.rpca"、 "harmony"、 "integrated.mnn" 、"umap.harmony" 、"tsne.harmony"'
),
  byrow=T,ncol=5
)


## 读取参数
args=getopt(command)

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name) || is.null(args$gmt)) {
  cat("Usage: Rscript cluster_umap -i input.rds -o outdir -n name\n")
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

if (!file.exists(args$gmt)) {
  cat("Error: gmt file not exists!\n")
  q()
}

if ( is.null(args$reduction)) {
  reduction <- "umap"
}else if (args$reduction == "umap"){
  reduction <- "umap"
}else if (args$reduction == "tsne"){
  reduction <- "tsne"
} else if (args$reduction == "integrated.cca"){
  reduction <- "integrated.cca"
} else if (args$reduction == "integrated.rpca"){
  reduction <- "integrated.rpca"
} else if (args$reduction == "harmony"){
  reduction <- "harmony"
} else if (args$reduction == "umap.harmony"){
  reduction <- "umap.harmony"
} else if (args$reduction == "tsne.harmony"){
  reduction <- "tsne.harmony"
}else if (args$reduction == "integrated.mnn"){
  reduction <- "integrated.mnn"
}else{
  cat("Error: reduction method not exists!\n")
  q()
}


library(Seurat)
library(qusage)
library(ggplot2)
library(scCustomize)
library(patchwork)
gene_set <- qusage::read.gmt(args$gmt)


print("读取RDS文件")
pbmc <- readRDS(args$input)
print("读取RDS文件完成")

sample_number <- length(unique(pbmc$orig.ident))
num_columns <- min(3, sample_number)


prefix <- paste0(args$outdir, "/by_gene")
dir.create(prefix, showWarnings = FALSE)

for (i in 1:length(gene_set)) {
  name =  names(gene_set[i])
  print( paste0("开始绘制", name, "的图"))
  for (j in 1:length(gene_set[[i]])) {
    if ( ! gene_set[[i]][j] %in% Features(pbmc) ){
      cat("Error: ", gene_set[[i]][j], "not in pbmc!\n")
    }else{
      print( paste0("开始绘制", name, "的", gene_set[[i]][j], "的图"))
      a_gene <- gene_set[[i]][j]
      height <- ceiling( sample_number / num_columns) * 4 

      print(paste0( "height: ", height))
      tryCatch({
        pdf(paste0(prefix, "/", name, "_" , a_gene, ".pdf") , width = 16, height = 18)
        p1 <- DimPlot_scCustom(pbmc, reduction = reduction)
        p2 <- FeaturePlot_scCustom(pbmc, features= a_gene , reduction = reduction) 
        p3<- VlnPlot_scCustom(pbmc, features = a_gene , raster = FALSE)
        p4 <- VlnPlot_scCustom(pbmc, features = a_gene , split.by = "orig.ident", raster = FALSE)
        

        p <- (p1 | p2) / p3 / p4 
        print(p)
        dev.off()

        pdf(paste0(prefix, "/", name, "_" , a_gene, "_by_sample.pdf") , width = 9, height = height)
        p2 <- FeaturePlot_scCustom(pbmc, features= a_gene , reduction = reduction, split.by = "orig.ident" , num_columns = num_columns, repel = TRUE)

        print(p2)
        dev.off()
      }, error = function(e) {
        print(e)
        cat("Error: ", gene_set[[i]][j], "not in pbmc!\n")
      })
    }
  }

}


### 绘制所有基因的DotPlot图
print("开始绘制所有基因的DotPlot图")
all_gene <- unique(Reduce(union, gene_set))
pdf(paste0(args$outdir, "/", args$name, "_", "all_gene_dotplot", ".pdf") , width = 10, height = 10)
p1 <- DotPlot_scCustom(pbmc, features = all_gene ,  flip_axes = T,   remove_axis_titles = FALSE) 
print(p1)
if ( length(all_gene) > 10){
	p2 <- Clustered_DotPlot(pbmc , features = all_gene)
	print(p2)
}
dev.off()

### 绘制所有基因的FeaturePlot图
print("开始绘制所有基因的FeaturePlot图")
pdf(paste0(args$outdir, "/", args$name, "_", "all_gene_featureplot", ".pdf") , width = 10, height = 10)
p1 <- DimPlot_scCustom(pbmc, reduction = reduction)
print(p1)

gene_number <- length(all_gene)
num_columns <- min(3, gene_number)
p2 <- FeaturePlot_scCustom(pbmc, features = all_gene , reduction = reduction,num_columns = num_columns, repel = TRUE)
print(p2)
dev.off()


### 绘制每个样本的UMAP图
all_samples <- unique(pbmc$orig.ident)
print("开始绘制所有样本的FeaturePlot图")
prefix <- paste0(args$outdir, "/by_sample") 
dir.create(prefix, showWarnings = FALSE)

for (i in all_samples) {
  print(paste0("开始绘制", i, "的FeaturePlot图"))
  pbmc_subset <- subset(pbmc, subset = orig.ident== i)
  prefix_of_sample <- paste0(prefix, "/", i)
  dir.create(prefix_of_sample, showWarnings = FALSE)

  print(paste0("开始绘制", i, "的UMAP图"))
  pdf(paste0(prefix_of_sample, "/", i, "_umap.pdf") , width = 8, height = 6)
  p1 <- DimPlot_scCustom(pbmc_subset, reduction = reduction) + ggtitle( paste0("UMAP of ",  i))
  print(p1)
  dev.off()
  
  print(paste0("开始绘制", i, "的FeaturePlot图"))
  height <- ceiling(length(all_gene) / 3) * 4 
  
  pdf(paste0(prefix_of_sample, "/", i, "_featureplot", ".pdf") , width = 12, height = height)
  p3 <- FeaturePlot_scCustom(pbmc_subset, features = all_gene , reduction = reduction , num_columns = num_columns, repel = TRUE) + ggtitle( paste0("FeaturePlot of ",  i))
  #p <- wrap_plots(p1, p3, ncol = 1)
  print(p3)
  dev.off()
}

### 绘制每个基因集的HeatMap
print("开始绘制每个基因集的HeatMap")
prefix <- paste0(args$outdir, "/by_gene_set")
dir.create(prefix, showWarnings = FALSE)
for (i in 1:length(gene_set)) {
  name =  names(gene_set[i])
  print( paste0("开始绘制", name, "的HeatMap"))
  pdf(paste0(prefix, "/", name, "_heatmap.pdf") , width = 10, height = 10)
  p1 <- DoHeatmap(pbmc, features = gene_set[[i]],  assay = args$assay, label = T, group.bar = T, raster = F)
  print(p1)

  dev.off()
}


 
