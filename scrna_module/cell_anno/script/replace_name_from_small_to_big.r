
library(getopt)
 
command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'small', 's', 1,'character', '小样本名',
  'big', 'b', 1,'character', '大样本名',
  'cluster', 'c', 1,'character', '小rds中cluster的列名， 默认为cluster',
  'subcell', 'l', 1,'character', '大rdss中subcell的列名， 默认为sub_celltype',  
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名',
  'misscellname', 'm', 1,'character', '没有对应上的subcelltype，默认为unknown, 设置usebig时，使用大rds中Idents 列名来填充， 否则使用unknown'
   ),

  byrow=T,ncol=5
  )

## 读取参数
args=getopt(command)

if (!is.null(args$help) || is.null(args$big) || is.null(args$outdir) || is.null(args$name) || is.null(args$small) ) {
  cat("Usage: Rscript this.r -s samll.rds -b big.rds -o outdir -n name\n")
  cat(paste(getopt(command, usage = T), "\n"))
  q(status = 1)
}


infile <- args$input
outdir <- args$outdir

if (!file.exists(outdir)) {
  dir.create(outdir)
}

name <- args$name

prefix <- paste(outdir, name, sep = "/")

small_cluster_name <- "cluster"
if (!is.null(args$cluster)) {
  small_cluster_name <- args$cluster
}
big_subcell_name <- "sub_celltype"
if (!is.null(args$subcell)) {
  big_subcell_name <- args$subcell
}

library(Seurat)
library(dplyr)

pbmc_small <- readRDS(args$small)
pbmc_big <- readRDS(args$big)

## check big_subcell_name 是否存在
if (big_subcell_name %in% colnames(pbmc_big@meta.data)) {
    cat(paste0("big_subcell_name:", big_subcell_name, " is in the big rds file , please check it again!"))
    cat("\n")
    print( colnames(pbmc_big@meta.data))
} 

if ( ! small_cluster_name %in% colnames(pbmc_small@meta.data)) {
    cat(paste0("small_cluster_name:", small_cluster_name, " is not in the small rds file , please check it again!"))
    cat("\n")
    print( colnames(pbmc_small@meta.data))
    stop()
}

## 首先给big 加上unknown

if (is.null(args$misscellname)) {
  pbmc_big@meta.data[[big_subcell_name]] <- "unknown"
}else if (args$misscellname == "usebig") {
  pbmc_big@meta.data[[big_subcell_name]] <- Idents(pbmc_big)
}else{
  cat("-m is not supported, please check it again! , use -m usebig")
  stop()
}


name_df <- data.frame(name = rownames(pbmc_small@meta.data) , cluster = pbmc_small@meta.data[[small_cluster_name]])

cell_indices <- match(name_df$name, colnames(pbmc_big))
pbmc_big@meta.data[[big_subcell_name]][cell_indices] <- name_df$cluster


print("save the result")
saveRDS(pbmc_big, paste0(prefix, ".rds"))