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
#read config.ini
ini<-opt$config
#ini.list <- read.config(file = ini)


outdir=opt$outdir #"../result/2_monocle3/"
covid_pb<-readRDS(opt$indir)
#构造cds对象
#基因表达矩阵-稀疏矩阵，expression_data和gene_anno的行是对应的，expression_data的列宇cell_metadata的行对应
expression_data <- GetAssayData( covid_pb, assay='RNA',slot='counts')
#数据框，行为细胞，列为细胞的属性
cell_metadata <- covid_pb@meta.data
#数据框，行为feature(基因)，列为基因属性
gene_anno <- data.frame(gene_short_name=rownames(expression_data))
rownames(gene_anno) <- rownames(expression_data)

#创建cds对象
cds <- new_cell_data_set(expression_data=expression_data, cell_metadata=cell_metadata, gene_metadata=gene_anno)


ini.list_Para_seurat_clusters="seurat_clusters"
ini.list_Para_seurat_title ="clusters"
ini.list_Para_cores=6
ini.list_Para_resolution=1e-3

prefix=opt$prefix
setwd(outdir)
###导入seurat对象
###monocle3
mkdirs(outdir,'1_trajectory/')
setwd(paste(outdir,'1_trajectory/',sep='/'))

if(ini.list_Para_seurat_clusters == 'seurat_clusters'){
	print("按照cluster进行分析，不需要输入-n参数")
}else{
	print("按照细胞类型进行分析，需要输入-n参数")
	if(is.null(opt$cellname)){
		print("Please give the celltype_marker.xls file ...")
		quit()
	}else{
		cell_name<-read.table(opt$cellname,header=T,sep='\t',comment.char='')
		cell_Group <- subset(cell_name,Group==opt$prefix)[,2:3]
		colnames(cell_Group) <- c('seurat_clusters', 'Cell_type')
		labers=cell_Group[match(as.character(covid_pb@meta.data$seurat_clusters),cell_Group[,1]),2]
		covid_pb$Cell_type <-labers
	}
}


###Monocle3聚类分区
print("预处理cds")
cds <- preprocess_cds( cds, method = "PCA",num_dim=100 )
cds <- reduce_dimension(cds, preprocess_method="PCA",reduction_method=c("UMAP"))
cds <- cluster_cells(cds)
#pdf(paste(prefix,"cell_type.fine.partition.pdf",sep='_'),w=8,h=4)
pdf(paste(prefix,"cell_type.fine.partition.pdf",sep='_'),w=20,h=8)

#p1 <- plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by =ini.list_Para_seurat_clusters)+ggtitle(ini.list_Para_seurat_title)
#p2 <- plot_cells(cds, color_cells_by = "partition", group_cells_by="partition",show_trajectory_graph = FALSE)+ggtitle("lable by partitionID")
p1 <- plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by =ini.list_Para_seurat_clusters,label_cell_groups = FALSE)+ggtitle(ini.list_Para_seurat_title)+theme(legend.position="right")
p2 <- plot_cells(cds, color_cells_by = "partition", group_cells_by="partition",show_trajectory_graph = FALSE,label_cell_groups = FALSE)+ggtitle("lable by partitionID")+theme(legend.position="right")
p3<-plot_grid(p1, p2)
print(p3)
dev.off()
##识别轨迹
cds <- learn_graph(cds)
cds2 <- learn_graph(cds,use_partition = F)
pdf(paste(prefix,'learn_graph_cells.pdf',sep='_'),w=20,h=8)
p1 <- plot_cells(cds, color_cells_by = 'cluster',label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)+theme(legend.position="right")
#p2 <- plot_cells(cds, color_cells_by = ini.list_Para_seurat_clusters, label_leaves = FALSE, label_branch_points = FALSE)+theme(legend.position="top")
p2 <- plot_cells(cds, color_cells_by = ini.list_Para_seurat_clusters, label_leaves = FALSE, label_branch_points = FALSE)+theme(legend.position="right") #去除图例

p3<-plot_grid(p1, p2)
print(p3)
dev.off()
pdf(paste(prefix,'learn_graph_cells.no.use_partition.pdf',sep='_'),w=20,h=8)
p1 <- plot_cells(cds2, color_cells_by = 'cluster',label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)+theme(legend.position="right")
#p2 <- plot_cells(cds2, color_cells_by = ini.list_Para_seurat_clusters, label_leaves = FALSE, label_branch_points = FALSE)+theme(legend.position="top")
p2 <- plot_cells(cds2, color_cells_by = ini.list_Para_seurat_clusters, label_leaves = FALSE, label_branch_points = FALSE)+theme(legend.position="right") #去除图例
p3 <- plot_cells(cds2, color_cells_by = 'orig.ident',label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)+theme(legend.position="right")
p4<-plot_grid(p1, p2, p3)
print(p4)
dev.off()

#mess=paste("细胞总数",format(nrow(tmp@meta.data)),"小于config.ini文件里面的maxcellnum的数字，不进行细胞数目的截取")
get_earliest_principal_node <- function(cds,colname,time_bin="Class-switched B"){
  cell_ids <- which(colData(cds)[, colname] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

if (is.null(opt$root)){
	cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds,ini.list_Para_seurat_clusters,time_bin=covid_pb@meta.data[,ini.list_Para_seurat_clusters][1]))
	cds2 <- order_cells(cds2, root_pr_nodes=get_earliest_principal_node(cds2,ini.list_Para_seurat_clusters,time_bin=covid_pb@meta.data[,ini.list_Para_seurat_clusters][1]))
}else{
	cds <- order_cells(cds, root_cells=opt$root)
	cds2 <- order_cells(cds2, root_cells=opt$root)
}
pdf(paste(prefix,'order_cells_cells.no.use_partition.pdf',sep='_'),w=8,h=4)
p1 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,label_branch_points = FALSE)
p2 <- plot_cells(cds2, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,label_branch_points = FALSE)
p3 <- plot_cells(cds, color_cells_by = "orig.ident", label_cell_groups = FALSE, label_leaves = FALSE,label_branch_points = FALSE)
p4 <- plot_cells(cds2, color_cells_by = "orig.ident", label_cell_groups = FALSE, label_leaves = FALSE,label_branch_points = FALSE)

p<-plot_grid(p1, p2, p3, p4)
print(p)
dev.off()


mkdirs(outdir,'2_track_gene/')
setwd(paste(outdir,'2_track_gene/',sep='/'))
##寻找拟时轨迹差异基因
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
#空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
cds <- estimate_size_factors(cds)
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=as.numeric(ini.list_Para_cores))
Track_genes['gene_short_name']=rownames(Track_genes)
###status p_value morans_test_statistic morans_I vst.mean vst.variance vst.variance.expected vst.variance.standardized vst.variable q_value gene_short_name
###status p_value morans_test_statistic morans_I gene_short_name q_value

Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.table(Track_genes,paste(outdir,'2_track_gene/',paste(prefix,"Trajectory_DEG_genes.xls",sep='_'),sep='/'),sep="\t",quote=F, row.names = F)
write.table(head(Track_genes),paste(outdir,'2_track_gene/',paste(prefix,"Trajectory_DEG_genes.example.xls",sep='_'),sep='/'),sep="\t",quote=F, row.names = F)
##差异基因表达趋势图
pdf(paste(outdir,'2_track_gene/',paste(prefix,'Genes_Jitterplot.pdf',sep='_'),sep="/"),width = 8, height = 6)
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by=ini.list_Para_seurat_clusters, 
                             min_expr=0.5, ncol = 2,label_by_short_name = FALSE)
print(p)
dev.off()
##绘制FeaturePlot图
#pdf(paste(outdir,paste(prefix,'Genes_Featureplot.pdf',sep='_'),sep="/"),width = 20, height = 6)
#p <- plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,label_cell_groups=FALSE,  label_leaves=FALSE)
#p$facet$params$ncol <- 5
#print(p)
#dev.off()
mkdirs(outdir,'3_cogene_model/')
setwd(paste(outdir,'3_cogene_model/',sep='/'))
##寻找共表达基因模块
#Track_genes <- read.csv("Trajectory_genes.csv")
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
cds<-preprocess_cds(cds, num_dim = 100)
cds@reduce_dim_aux$gene_loadings <- covid_pb@reductions[["pca"]]@feature.loadings
if(length(genelist) < 2){
	print("graph_test() 函数最终得到的基因数量小于2，做不了共表达模块的分析")
	q()
}else{
	print("基因数量大于2 ，可以做共表达模块的分析")
}
gene_module <- find_gene_modules(cds[genelist,], resolution=as.numeric(ini.list_Para_resolution), cores = as.numeric(ini.list_Para_cores))
write.table(gene_module,paste(outdir,'3_cogene_model/',paste(prefix,"Genes_Module.xls",sep='_'),sep='/'),sep="\t",quote=F, row.names = F)

pdf(paste(outdir,'3_cogene_model/',paste(prefix,'Genes_Module.pdf',sep='_'),sep="/"),width = 8, height = 8)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)[,ini.list_Para_seurat_clusters])
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module", row.names(agg_mat))
p <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
print(p)
dev.off()
saveRDS(cds, file = paste(outdir,'3_cogene_model/',paste(prefix,'immune_combined_monocle3.rds',sep='_'),sep="/"))



