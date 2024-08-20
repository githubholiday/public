
library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名',
  "group_by" , 'c', 2, 'character', '聚类的名称，默认为"seurat_clusters"',
  "split_by" , 's', 2, 'character', '按照哪个列进行拆分，默认为orig.ident'
),
  byrow=T,ncol=5
)

# 读取参数
args=getopt(command)

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name)) {
  cat("Usage: Rscript cluster_umap -i input.rds -o outdir -n name\n")
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
if ( is.null(args$group_by)) {
  cluster_name <- "seurat_clusters"
}else{
    cluster_name <- args$cluster_name
}

if ( is.null(args$split_by)) {
  split_by <- "orig.ident"
}else{
    split_by <- args$split_by
}




library(SCpubr)
library(Seurat)


pbmc <- readRDS(args$input)

prefix <- paste0(args$outdir, "/", args$name)

pdf(paste0(prefix, "_cluster_composition.pdf") , width = 8, height = 6)
p1 <- SCpubr::do_BarPlot(sample = pbmc, 
                         group.by = cluster_name, 
                         legend.position = "none", 
                         plot.title = "Number of cells per cluster")
print(p1)
dev.off()

pdf(paste0(prefix, "_cluster_composition_by_sample.pdf") , width = 8, height = 6)
p1 <- SCpubr::do_BarPlot(sample = pbmc,
                         group.by = cluster_name,
                         split.by = split_by, 
                         legend.position = "right",
                         plot.title = "Number of cells per cluster in each sample",
                         position = "stack")

p2 <- SCpubr::do_BarPlot(sample = pbmc,
                         group.by = split_by,
                         split.by = cluster_name,
                         legend.position = "right",
                         plot.title = "Number of cells per sample in each cluster",
                         position = "stack")

print(p1)
print(p2)

dev.off()

pdf(paste0(prefix, "_cluster_composition_by_sample_percentage.pdf") , width = 8, height = 6)
p1 <- SCpubr::do_BarPlot(sample = pbmc,
                         group.by = cluster_name,
                         split.by = split_by, 
                         legend.position = "right",
                         plot.title = "Pecentage of cells per cluster in each sample",
                         position = "fill")
print(p1)


p2 <- SCpubr::do_BarPlot(sample = pbmc,
                         group.by = split_by,
                         split.by = cluster_name,
                         legend.position = "right",
                         plot.title = "Pecentage of cells per sample in each cluster",
                         position = "fill")

print(p2)
dev.off()

output_xls <- paste0(prefix, "_cluster_composition.xls")
new_df <- data.frame(
  cluster_name = pbmc[[cluster_name]],
  sample_name = pbmc[[split_by]]

)
result <- table( new_df )
write.table(result , output_xls , sep = "\t", quote = F)
print("细胞比例统计结果已保存")

