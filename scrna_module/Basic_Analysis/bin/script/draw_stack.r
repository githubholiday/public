library(getopt)

command=matrix(c( 
    'help', 'h', 0,'logical', '帮助文档',
    'rds', 'r', 1, 'character', '输入的rds文件',
    'outdir', 'o', 1, 'character', '输出的目录',
    'prefix' , 'p', 1, 'character', '输出的文件名',
    'Xname' , 'X', 1, 'character', '用于绘图的X轴的slots名称,一般是orig.ident,group,Phase',
    'Yname' , 'Y', 1, 'character', '用于绘图的Y轴的slots名称,一般是seurat_clusters,CellType'
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
    -X:重命名前的slot名称
    -Y:重命名后的slot名称

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

annotation_info <- Pull_Cluster_Annotation(annotation = args$input)
pbmc <-  Rename_Clusters(seurat_object = pbmc_t, new_idents = annotation_info$new_cluster_idents)

#如果更换名称的话，可以改CellType
pbmc[[new_idents]] <- Idents(pbmc)

rds_file <- paste(out_pre, ".rename.rds", sep="")
saveRDS( pbmc, rds_file)

#绘图
pdf_file <- paste(out_pre, ".dimplot.pdf", sep="")
pdf(pdf_file)
#a<- DiscretePalette_scCustomize(num_colors = 26, palette = "alphabet")
p1 <- DimPlot_scCustom(pbmc, reduction = "umap.harmony", label = TRUE)# , colors_use =a)
#p2 <- DimPlot_scCustom(pbmc_new, reduction = "tsne", label = TRUE , colors_use = a)

p2 <- DimPlot_scCustom(pbmc, reduction = "umap.harmony", label = TRUE , group.by = "seurat_clusters" ) #, colors_use = a)
#p4 <- DimPlot_scCustom(pbmc_new, reduction = "tsne", label = TRUE , group.by = "seurat_clusters" , colors_use = a)
print(p1)
print(p2)
dev.off()

# 输出表格以及如果样本多的话，绘制样本和组别的stackplot图
celltype_num <- as.data.frame(table(pbmc[[new_idents]]))
colnames(celltype_num)<-c("clusters","count")
count_file <- paste(out_pre, ".celltype_count.xls",sep="")
write.table( celltype_num ,count_file , sep="\t",quote=FALSE,row.names=FALSE)

# 堆叠图
tt <- table(pbmc@meta.data$orig.ident,pbmc[[new_idents]])
tt_new <- cbind(Sample=rownames(tt),tt)
write.table(tt_new,"HnEC/02_cellcycle/HnEC.cellcycle.xls",sep="\t",quote=FALSE, row.names=FALSE)
#计算细胞比例
Cellratio <- prop.table(table(pbmc@meta.data$Phase,pbmc@meta.data$CellType))
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("CellCycle","Cluster","Freq")

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
library(ggplot2)
pdf("HnEC/02_cellcycle/HnEC_cellcycle_stackplot.pdf")
p1<-ggplot(Cellratio) + 
  geom_bar(aes(x =Cluster, y= Freq, fill = CellCycle),stat = "identity",position="fill",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Cluster',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
print(p1)
dev.off()


