library('getopt')
para<-matrix(c(
	'help',	'h',	0,	"logical",
	'datafile',	'd',	1,	"character",
	'outdir',	'o',	1,	"character",
	'sample',	's',	1,	"character",
	'testmethod',	'tm',	2,	"character",
	'Ident',	"I",2,"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=FALSE)

print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	Usage example:
	Rscript seurat_analysis.R -d rdsfile -s sample -o outdir/
	Options:
	--help	h	NULL	get this help
	--datafile	d	character	the rds file [forced]
	--sample	s	character	SampleName[forced]
	--outdir	o	character	output file dir [forced]
	--testmethod	tm	character	test method for findMakers [default: wilcox]
	\n")
	q(status=1)
}
if ( is.null(opt$datafile) || is.null(opt$sample) || is.null(opt$outdir)){ print_usage } 
if ( is.null(opt$testmethod))	{ opt$testmethod <- c("wilcox") }

.libPaths()
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(reticulate)
library(SeuratWrappers)
library(Cairo)
library(gridExtra)

prefix <- paste(opt$outdir,opt$sample,sep='/')


## creat object
object <- readRDS(opt$datafile)
object@assays
#Idents(object) <- opt$Ident

# Expression Violin plot
Spot_exp_plot <- paste(prefix,'Spots_Count.pdf',sep='_')
pdf(Spot_exp_plot,width=16,height=8) ## 注意画图输出
plot1 <- VlnPlot(object = object, features = "nFeature_Spatial", ncol = 1, combine = "False")
plot2 <- VlnPlot(object = object, features = "nCount_Spatial", ncol = 1, combine = "False")
CombinePlots(plots=c(plot1,plot2), ncol=2) #输出空间表达小提琴图
print(Spot_exp_plot)
dev.off()


# UMI correlation scatter plot
graph<-paste(prefix,'CorPlot.pdf',sep='_')
pdf(graph,width=12,height=8)
plot1 <- FeatureScatter(object = object, feature1 = "nFeature_Spatial", feature2 = "nCount_Spatial")
plot1 #输出feature和count相关性
print(graph)
dev.off()


## 空间Spots覆盖深度显示
graph<-paste(prefix,'nCount_Spatial.pdf',sep='_')
pdf(graph,width=16,height=8)
plot1 <- VlnPlot(object, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(object, features = "nCount_Spatial") + theme(legend.position = "right")
plot3 <- plot_grid(plot1, plot2) ##输出空间表达图
plot3
print(graph)
dev.off()
ggsave(plot3, filename=paste(prefix,'nCount_Spatial.png',sep='_'), width = 16, height = 8)


## 写矩阵文件
file<-paste(prefix,'all_spots.csv',sep='_')
print(file)
write.csv(object@meta.data,file=file,quote=F)

#file<-paste(prefix,'all_UMI.csv',sep='_')
#print(file)
#:wwrite.csv(object@assays$Spatial@data,file=file,quote=F)


## SCT标准化
object <- subset(object, nCount_Spatial>0)
object <- SCTransform(object, assay = "Spatial", verbose = FALSE)


## PCA heatmap 
cell_for_heatmap = round(length(object@meta.data[,1])*0.05)
data <- object@assays$SCT@scale.data
file <- paste(prefix,'Norm_Scale.csv',sep='_')
print(file)
write.csv(data,file=file,quote=F)


## Find VariableGenes for PCA
graph<- paste(prefix,'VariableGenes.pdf',sep='_')
pdf(graph,w=12,h=8)

#object <- FindVariableFeatures(object = object, mean.function = "FastExpMean", dispersion.function = "FastLogVMR", x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff= 0.5, do.contour = FALSE) ## calculates highly variable genes and focuses on these for downstream analysis，为PCA做准备  ==== add do.contour = FALSE 20180926 by yaomengcheng 只会对画图有影响， Draw contour lines calculated based on all genes
top10 <- head(VariableFeatures(object), 10)
plot1 <- VariableFeaturePlot(object)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(graph)
dev.off()


data <- object@assays$SCT@var.features 
file <- paste(prefix,'VariableGenes.csv',sep='_')
print(file)
write.csv(data,file=file,quote=F)


#object <- ScaleData(object = object, vars.to.regress="percent.mt")
# dimensional
object <- RunPCA(object = object, assay = "SCT", verbose = FALSE, features = VariableFeatures(object = object), npcs=20)
object <- FindNeighbors(object, reduction = "pca", dims = 1:20)

#对不同分辨率的cluster个数<20进行动态调整resolution
values <- seq(from = 0.8, to = 0.2, by = -0.2)
for (value in values) {
	object1 <- FindClusters(object, verbose = FALSE, resolution = value)
	cluster_count <- length(levels(object1$seurat_clusters))
	print(paste("resolution:",value,"cluster count",cluster_count,sep=" "))
	if (cluster_count < 12){
		print(paste("采用此分辨率进行后续分析",value, sep=":"))
	}
    break
}

celltype_num <- as.data.frame(table(object@meta.data$seurat_clusters))
colnames(celltype_num)<-c("clusters","count")
outfile <- paste(prefix , "_celltype_count.xls",sep="")
write.table( celltype_num ,  outfile, sep="\t",quote=FALSE,row.names=FALSE)

#object <- FindClusters(object, verbose = FALSE, resolution = 0.5)
object <- RunUMAP(object1, dims = 1:20, metric="correlation")


### plot cluster 
graph <- paste(prefix,'cluster_umap.pdf',sep='_')
print(graph) 
pdf(graph,w=16,h=8)
p1 <- DimPlot(object, reduction = "umap", label = TRUE, raster = FALSE)
p2 <- SpatialDimPlot(object, label = TRUE, label.size = 1, stroke = NA)
p3 <- plot_grid(p1, p2)
p3
dev.off()
ggsave(p3, filename=paste(prefix,'cluster_umap.png',sep='_'), width = 16, height = 8)


#PCA plot
graph<- paste(prefix,'PCA.pdf',sep='_')
pdf(graph,w=12,h=8)
DimPlot(object = object, reduction="pca")
print(graph)
dev.off()


### plot violn for cluster
graph <- paste(prefix,'Vln_cluster.pdf',sep='_')
pdf(graph,width=14,height=10)
p1<-VlnPlot(object = object, features = "nFeature_Spatial", ncol = 1, combine = "False", pt.size = 0, group.by="seurat_clusters")
p2<-VlnPlot(object = object, features = "nCount_Spatial", ncol = 1, combine = "False", pt.size = 0, group.by="seurat_clusters")
CombinePlots(plots=c(p1,p2), ncol=1)
print(graph)
dev.off()


##find markers for every cluster compared to all remaining cells, report only the positive ones 
object.markers <- RunPrestoAll(object,only.pos = FALSE, min.pct =0.20, logfc.threshold = 0.2,test.use=opt$testmethod)

#object.markers <- FindAllMarkers(object = object, only.pos = FALSE, min.pct = 0.20, logfc.threshold = 0.25, test.use =opt$testmethod)
AllMakers <- paste(prefix,'all_markers_result.csv',sep='_')
marker_result <- subset(object.markers, select=c(8,1,3,4,5,6,7))
#marker_result <- subset(object.markers, select=c(7,1,2,3,4,5,6))
names(marker_result)[names(marker_result)=='gene'] <- 'gene_name'
write.csv(marker_result,file=AllMakers,quote=F,row.names=F)

#输出每个cluster的marker基因
for (i in unique(Idents(object))){
	markers <- marker_result[ marker_result$cluster==i,]
	outfile <- paste(prefix, i, "marker.xls",sep="_")
	write.table(markers, outfile, sep="\t", quote=FALSE,row.names = FALSE )
    }


object@misc$markers <- object.markers

file <- paste(prefix,'all_analysis.rds',sep='_')
saveRDS(object, file)

