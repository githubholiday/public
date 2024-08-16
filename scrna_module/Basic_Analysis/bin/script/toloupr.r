library(getopt)
command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名'
  ),
  byrow=T,ncol=5
)

## 读取参数
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


library(Seurat)
library(loupeR)


cat(" 如果是第一次跑，请运行 library(loupeR);loupeR::setup()") 

output_name <- paste0(args$outdir, "/", args$name, ".cloupe")
pbmc <- readRDS(args$input)
create_loupe_from_seurat(pbmc , output_name= output_name, force = TRUE)
print("Done!")



