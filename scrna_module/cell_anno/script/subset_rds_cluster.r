library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入subset后的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名' ,
  'resolution' , 'r', 1, 'numeric', '分辨率 默认为 0.6'
),
  byrow=T,ncol=5
)

## 读取参数
args=getopt(command)

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name) ) {
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

resolution <- 0.6
if (!is.null(args$resolution)) {
  resolution <- args$resolution
}
print("聚类分辨率为"  )
print(resolution)


output_umap_result <- function(pbmc , prefix , method = "umap" , integrate=TRUE , seed = 0  ) {
    if (integrate){
        prefix <- paste0(prefix, "_integrate")
    }else {
        prefix <- paste0(prefix, "_unintegrate")
    }
    if (seed > 0){
        prefix <- paste0(prefix, "_seed_", seed)
    }

    pdf(file = paste0(prefix, "_" , method , ".pdf"))

    plot_name= paste0(method , "_plot")
    plot1 <- DimPlot(pbmc, reduction = method) + ggtitle(  paste0( plot_name ,  "of all samples"))
    print(plot1)
    plot11 <- DimPlot(pbmc, reduction = method , label=TRUE) + ggtitle( paste0( plot_name ,  "of all samples"))
    print(plot11)
    
    plot11 <- DimPlot(pbmc, reduction = method , split.by= "orig.ident") + ggtitle( paste0( plot_name ,  "of all samples"))
    print(plot11)
    plot2 <- DimPlot(pbmc, reduction = method, group.by="orig.ident") + ggtitle(paste0( plot_name ,  "of all samples"))
    print(plot2)
    if ("group" %in% colnames(pbmc)) { 
      plot2 <- DimPlot(pbmc, reduction = method, group.by="group") + ggtitle(paste0( plot_name ,  "of Group"))
      print(plot2)
      plot2 <- DimPlot(pbmc, reduction = method, split.by="group") + ggtitle(paste0( plot_name ,  "of Group"))
      print(plot2)
    }
    dev.off()
}

library(Seurat)
library(ggplot2)

pbmc <- readRDS(args$input)

pbmc<- FindNeighbors(pbmc, reduction = "harmony", dims = 1:30)

pbmc <- FindClusters(pbmc, resolution = resolution , cluster.name = "harmony_clusters" )

pbmc <- RunUMAP(pbmc,  reduction="harmony", reduction.name = "umap.harmony"  , dims = 1:30)

prefix <- paste0(args$outdir, "/", args$name , "_" , resolution )
output_umap_result(pbmc , prefix , method = "umap.harmony")

saveRDS(pbmc, paste0(prefix, ".umap.rds"))  

celltype_num <- as.data.frame(table(pbmc@meta.data$CellType))
colnames(celltype_num)<-c("CellType","count")
count_file <- paste(prefix, ".celltype_count.xls",sep="")
write.table( celltype_num ,count_file , sep="\t",quote=FALSE,row.names=FALSE)

print("Done")

# subset <- FindNeighbors(subset, reduction = "harmony", dims = 1:30)
# subset <- FindClusters(subset, resolution = 0.5 , cluster.name = "harmony_clusters" )
# subset <- RunUMAP(subset, reduction = "harmony", dims = 1:30 ,  reduction.name = "umap.harmony" )
# DimPlot(subset, reduction = "umap.harmony" , label = TRUE)



