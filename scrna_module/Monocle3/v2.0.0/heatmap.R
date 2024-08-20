#!/annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/Rscript
#名称：monocle3.r
#作者：姚盟成
#邮箱：mengchengyao@genome.cn
#时间：2020-12-3
#版本：v0.0.1
#用途：利用monocle3进行10x 数据进行发育轨迹分析，此程序暂时只支持使用seuratV3.0以上版本的输出rds文件进行分析
###说明：
#程序开发环境/annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/Rscript，需要指定R，指定包的路径,使用为seurat3.0版本以上
#===========================================================
library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'prefix',	'p',	1,	"character",
	'indir',	'i',	1,	"character",
	'root',	'r',	2,	"character",
	'outdir',	'o',	1,	"character",
	'cellname',	'n',	2,	"character",
	'groupnum',	'g',	2,	"character",
	'config',	'c',	1,	"character"
),byrow=TRUE,ncol=4)
#===========================================================
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	========================================================================================================================================
	indir输入seurat标准输出的rds文件:
	目前只接受seurat3.0以上版本的输出rds文件
	========================================================================================================================================
	prefix:the output prefix of files,such as pictures and excel.
	========================================================================================================================================
	root 指定rootcell
	========================================================================================================================================
	========================================================================================================================================
	outdir:outdir  of outputs,we will setwd(opt$outdir)
	Usage example:
	Rscript this.r -i rds -o outdir -p prefix
	Options:
	--help		h	NULL		get this help
	--indir	i	character	indir for expression file[forced]
	--root	r	character	config.ini file for group and other Para[forced]
	--outdir	o	character	The	resurt of out dir for analysis [forced]
	--prefix	p	character	the prefix for outputfiles [forced]
	--config	c	character	the config.ini file [forced]
	--groupnum	g	character	the *group.celltype.xls file [optional]
	--cellname	n	character	the celltype_marker.xls file [optional]
	\n")
	q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )	{ print_usage(para) }
if ( is.null(opt$indir) )	{ cat("Please input the data file1 ...\n\n") ; print_usage(para)}
#if ( is.null(opt$root) )	{opt$database<-""}
if ( is.null(opt$outdir) )	{ cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$prefix) )	{ cat("Please give the prefix for outputfiles ...\n\n") ; print_usage(para) }
if ( is.null(opt$config) )	{ cat("Please give the config.ini file ...\n\n") ; print_usage(para) }

##这个分析用最新的seurat包进行分析

mkdirs <- function(outdir,fp) {
	if(!file.exists(file.path(outdir,fp))) {
#		mkdirs(dirname(fp))
		dir.create(file.path(outdir,fp))
	}else{
			print(paste(fp,"Dir already exists!",sep="     "))
		}
}

library(Seurat)
library(sctransform)
library(dplyr)
library(monocle3)
#library(monocle)
library(ggplot2)
library(cowplot)
#library(configr)
outdir <- opt$outdir
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
function_R <- paste(script.basename, 'function.R', sep='/')
source(function_R)

mkdirs <- function(outdir,fp) {
	if(!file.exists(file.path(outdir,fp))) {
		dir.create(file.path(outdir,fp))
	}else{
		print(paste(fp,"Dir already exists!",sep="     "))
	}
}

##############################################
library(SeuratWrappers)
cds <- readRDS(opt$indir)
cds_subset <- as.cell_data_set(cds)
heatmap_dir <- paste(outdir, "4_heatmap",sep="/")
mkdirs(outdir, "4_heatmap")
pbmc$psd <- floor(cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]])

pdf(paste(heatmap_dir,"/",prefix,"_umap_dimplot.pdf",sep=""))
p1<-DimPlot(pbmc , group.by="psd" , reduction="umap.harmony")
p2<-DimPlot(pbmc , group.by="orig.ident" , reduction="umap.harmony")
p<-plot_grid(p1,p2)
print(p)
dev.off() 


## 寻找随伪时间变化的基因
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=2)
## 按贡献值排序，输出结果
library(dplyr)
pr_graph_test_res <- pr_graph_test_res %>% arrange(desc(abs(morans_I)))

choose_gene <- rownames(subset(pr_graph_test_res, q_value < 0.05))[1:10]

pdf(paste(heatmap_dir,"/",prefix,"_trajectory.pdf",sep=""))
p<-plot_cells(cds_subset, 
       genes=choose_gene, show_trajectory_graph=T, 
       label_cell_groups=FALSE,  
       label_leaves=FALSE,
       label_branch_points=FALSE
)
print(p)
dev.off()

##### 画热图，https://github.com/cole-trapnell-lab/monocle-release/issues/295
modulated_genes <- graph_test(cds_subset, neighbor_graph = "principal_graph", cores = 2)
genes <- row.names(subset(pr_graph_test_res, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(cds_subset)[match(genes,rownames(rowData(cds_subset))),order(pseudotime(cds_subset))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

pdf(paste(heatmap_dir,"/",prefix,"_heatmap.pdf",sep=""))
htkm <- draw(htkm)
htkc <- draw(hthc)
p <- plot_grid(htkm, htkc)
print(p)
dev.off()

clusters <- row_order(htkm)
for (i in 1:length(clusters)) {
  genes <- rownames(pt.matrix)[clusters[[i]]]
  cat("Cluster:", i, "Genes:", genes, "\n")
}


 for (i in 1:length(row_order(htkm))){
  if (i == 1) {
    clu <- t(t(row.names(pt.matrix[row_order(htkm)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")
 } else {
    clu <- t(t(row.names(pt.matrix[row_order(htkm)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}
write.table(out, file= "gene_clusters.txt", sep="\t", quote=F, row.names=FALSE)



