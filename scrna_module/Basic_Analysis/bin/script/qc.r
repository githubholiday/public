library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'mt' , 'm', 2, 'character', '线粒体基因的正则表达式 , 默认为^MT-',
  'cutoff' , 'c', 2, 'numeric', '过滤的上下限,线粒体cutoff,默认为c(200,4000, 20)',
  'input', 'i', 1, 'character', '输入的合并后的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名'),
  byrow=T,ncol=5
  )
  

## 读取参数
args=getopt(command)


if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name)) {
  cat("Usage: Rscript qc.r -i input.rds -o outdir -n name\n")
  cat(paste(getopt(command, usage = T), "\n"))
    q() 
}

if (!file.exists(args$input)) {
  cat("Error: input file not exists!\n")
  q()
}

if (!dir.exists(args$outdir)) {
  cat("Error: outdir not exists!\n")
  q()
}

if ( is.null(args$mt)) {
  args$mt <- "^MT-"
}

if (is.null(args$cutoff)) {
  cutoff <- c(200,4000, 20)
}else{
   cutoff <- as.numeric(strsplit(args$cutoff, ",")[[1]])
}
print("过滤阈值为"  )
print(cutoff)

pbmc <- readRDS( args$input )
library(dplyr)
library(Seurat)
library(patchwork)

### 对线粒体基因进行质控
pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = args$mt)

pdf(paste(args$outdir, paste (args$name, "before_qc.pdf" , sep="."), sep="/"))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

before_qc_table = as.data.frame.array(table(Idents(pbmc)))

### 对双细胞进行质控
pbmc <- subset(pbmc, subset = nFeature_RNA > cutoff[1] & nFeature_RNA < cutoff[2] & percent.mt < cutoff[3])

after_qc_table = as.data.frame.array(table(Idents(pbmc)))
colnames(before_qc_table) = c("before_qc")
colnames(after_qc_table) = c("after_qc")

pdf(paste(args$outdir, paste (args$name, "after_qc.pdf" , sep="."), sep="/"))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

saveRDS(pbmc, paste(args$outdir, paste (args$name, "after_qc.rds" , sep="."), sep="/"))

sample_names <- rownames(before_qc_table)

merge2table = data.frame( sample_names,  before_qc_table, after_qc_table)
write.table(merge2table, paste(args$outdir, paste (args$name, "qc_table.txt" , sep="."), sep="/"), sep="\t", quote=F, row.names=F, col.names=T)



