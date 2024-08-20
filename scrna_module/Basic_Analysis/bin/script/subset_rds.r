library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名' ,
  'cluster' , 'c', 1, 'character', '指定某个聚类，比如1，多个用逗号分隔',
  'Ident' , "I" , 2 , 'character' , '聚类的标签，默认为harmony_clusters， 可选为任意列名'
),
  byrow=T,ncol=5
)

## 读取参数
args=getopt(command)

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name) || is.null(args$cluster)) {
  cat("Usage: Rscript Pseudobulk.r -i input.rds -o outdir -n name \n")
  cat(paste(getopt(command, usage = T), "\n"))
    q() 
}

if (!file.exists(args$input)) {
  cat("Error: input file not exists!\n")
  q()
}

if (!dir.exists(args$outdir)) {
  dir.create(args$outdir , recursive = T)
}

choose_clusters <- strsplit(args$cluster, ",")[[1]]



library(Seurat)

pbmc <- readRDS(args$input)


if ( !is.null(args$Ident) ){
  Idents(pbmc) <-args$Ident
  print("使用指定的标签")
  print(args$Ident)
}

print(table(Idents(pbmc)))

my_indent <- args$Ident
subset <- subset(pbmc, idents = choose_clusters)
#levels(subset@meta.data$my_indent)<- choose_clusters
saveRDS(subset, paste(args$outdir, paste (args$name, "subset.rds" , sep="."), sep="/"))

print("Done")

# subset <- FindNeighbors(subset, reduction = "harmony", dims = 1:30)
# subset <- FindClusters(subset, resolution = 0.5 , cluster.name = "harmony_clusters" )
# subset <- RunUMAP(subset, reduction = "harmony", dims = 1:30 ,  reduction.name = "umap.harmony" )
# DimPlot(subset, reduction = "umap.harmony" , label = TRUE)



