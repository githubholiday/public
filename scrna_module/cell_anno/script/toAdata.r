library(getopt)
command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名',
  'count', 'c', 0, 'logical', '是否导出count，默认是normalized count'
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
  cat("Error: outdir not exists!\n")
  q()
}

bool_count <- FALSE
if ( !is.null(args$count )){
  bool_count <- TRUE
}



library(reticulate)
#library(scater)


# 方法二
library(Seurat)
library(sceasy)

pbmc <- readRDS(args$input)
outfile <- paste0(args$outdir,"/",args$name,".h5ad")
#pbmc <- readRDS("/work/share/acuhtwkcu9/tuchengfang/04_Project/20240105DXB01ZD002/03_Delete_Sample/result/cluster2/rat.rds")
pbmc<- JoinLayers(pbmc)

pbmc[["RNA"]] <- as(pbmc[["RNA"]], "Assay")
use_python("/software/conda/bin/python3")

if (bool_count){
  pbmc[["RNA"]]$data <- pbmc[["RNA"]]$counts
}

sceasy::convertFormat(pbmc, from="seurat", to="anndata",  outFile=outfile , drop_single_values=FALSE)

print('Done!')


#Please install Python with one of following methods:
#- https://github.com/rstudio/python-builds/
#- reticulate::install_python(version = '<version>')




## 方法1 
#library(Seurat)
#pbmc<- JoinLayers(pbmc)
#pbmc.sce <-  as.SingleCellExperiment(pbmc)
#BiocManager::install("zellkonverter")
#library(zellkonverter)
#adata <- SCE2AnnData(pbmc.sce)
#out_path <- tempfile(pattern = "test.h5ad")
#writeH5AD(sce_zeisel, file = out_path)


#sceasy::convertFormat(pbmc.sce, from="sce", to="anndata",      outFile='filename.h5ad')




## 方法三 用 SeuratDisk 报错
#library(Seurat)
#library(SeuratData)
#library(SeuratDisk)
# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html



#rds$Sample = rds$stim
#DefaultAssay(rds)<-"RNA"
#Idents(rds) = rds$celltype
#rds$celltype=as.character(rds$celltype)
#rds@assays$RNA@scale.data = as.matrix(rds@assays$RNA@counts)
#SaveH5Seurat(rds, filename = "file.h5Seurat",overwrite=T)
#Convert("file.h5Seurat", dest = "h5ad")


#SaveH5Seurat(pbmc, filename = "pbmc.h5Seurat")
#Convert("pbmc.h5Seurat", dest = "h5ad")

##gplots::balloonplot(table(sce.all$orig.ident,Idents(pbmc)))
