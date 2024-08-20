library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名',
  'method' , 's', 2, 'character', '合并方法，默认是sum，可选mean',
  'colname' , 'c', 2, 'character', '列名，默认是orig.ident，可选为任意列名，也可以是多列，用逗号分隔',
  'assay', 'a', 2, 'character', 'Assay列，默认为RNA,可选Spatial'
  ),
  byrow=T,ncol=5
)

## 读取参数
args=getopt(command)

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name)) {
  cat("Usage: Rscript Pseudobulk.r -i input.rds -o outdir -n name\n")
  cat(paste(getopt(command, usage = T), "\n"))
    q() 
}

if (!file.exists(args$input)) {
  cat("Error: input file not exists!\n")
  q()
}

if (!dir.exists(args$outdir)) {
  cat("Error: outdir not exists!\n")
  dir.create(args$outdir , recursive = T)

}

if (is.null(args$assay)){
	args$assay <- "RNA"
}

if ( is.null(args$method)) {
  method <- "sum"
}else if (args$method == "sum"){
    method <- "sum"
}else if (args$method == "mean"){
    method <- "mean"
}else{
    cat("Error: method not exists!\n")
    q()
}

if ( is.null(args$colname)) {
  colname <- c("orig.ident")
}else{
    colname <- strsplit(args$colname, ",")[[1]]
}

library(Seurat)
library(Matrix.utils)
#library(tidyverse)

print("读取RDS"  )
pbmc <- readRDS(args$input)
## 合并多个layers
print("合并多个layers"  )
merged <- pbmc

tryCatch({
  merged <- JoinLayers(pbmc)
},
error = function(e) {
  print("没有多个layers"  )
})

print( paste0("使用" , colname))
if (length(colname) == 1){
  groups <- merged[[colname]]
}else{
  allgroups <- merged[[colname]]
  groups <- merged[[colname]][,1]
  for (i in 2:length(colname)){
    groups <- paste(groups, merged[[colname]][,i], sep = "_")
  }
}

print(head(groups))


counts <- merged[[args$assay]]$counts

aggr_counts <- aggregate.Matrix(t(counts),  groupings = groups, fun = method) 
aggr_counts <- t(aggr_counts)

## 生成矩阵
aggr_m <- as.matrix( aggr_counts )
outfile <- paste0(args$outdir, "/", args$name, ".xls")
write.table(aggr_m, file = outfile, sep="\t", quote=FALSE,  row.names = TRUE)


## grr
## install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)
## install.packages("tidyverse")

#matrix2<- read.table(paste0(args$outdir, "/", args$name, ".csv"))
#mat_sparse <- as(matrix2, "dgCMatrix")

