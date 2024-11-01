
library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名',
  'gmt' , 'g', 1, 'character', 'gmt文件',
  'assay', 'a', 1, 'character', 'assay 标签名，可选RNA Spatial SCT',
  'drawbysample' , 's', 2, 'character', '绘制每个样本的UMAP图，默认为FALSE，不绘制，可选为"T"、"F"',
  'reduction' , 'r', 2, 'character', '降维的方法，默认为"umap", 可选为"umap"、"tsne" 、"integrated.cca"、"integrated.rpca"、 "harmony"、 "integrated.mnn" 、"umap.harmony" 、"tsne.harmony"'
),
  byrow=T,ncol=5
)


## 读取参数
args=getopt(command)

print_usage <- function(para=NULL){
    cat("参数说明：
    -i:rds文件
    -o:输出目录，可以不存在
    -n:输出前缀
    -g:gene列表文件，以.gmt结尾 至少三列：geneset_name gene_set_description gene1 gene2
    -a:数据类型，assay的标签名，可选RNA，Spatial
    -r:降为方法，在画图的时候使用默认为umap,可选为map,tsne,integrated.cca,integrated.rpca,harmony,integrated.mnn,umap.harmony,tsne.harmony

    使用实例：
    /public/software/apps/singularity/3.7.3/bin/singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript draw_plot.r -i rds -o outdir -n sample -g gene.gmt -a Spatial 


    结果文件说明：
    1. by_gene
    1.1 by_gene/gene_.pdf：umap图，featureplot图，在不同群里的violinplot（comine)
    1.2 by_gene/gene_by_sample.pdf：在不同样本中的featureplot和violinplot
    2. by_gene_set/gene_set_dotplot.pdf:某个gene_set中所有基因在combine样本不同clsuter上的dotplot图
        3. by_sample：目前没画过，后续补充
    4. spatial
    4.1 featureplt_spatial.pdf: 基因在空间上的featureplot图
    4.2 spatial.pdf：是空间上的Umap和空间dimplot

    5. all_gene_dotplot.pdf：所有基因的dotplot图，可能会由于基因多，导致比较堆叠，可以删除
    6. all_gene_featureplot.pdf：所有基因的featureplot图，可能会由于基因多，导致比较堆叠，可以删除

    升级记录
    v3:2024-9-11 tx
    1. 如果是单个样本的不再绘制by_sample类型的图
    2. 增加by_gene_set的dotplot图 
    3. 增加spatial目录下按照gene_set绘制的空间featureplo图
    
    v4:2024-09-29 tx
    1. 将spatial目录下按照gene_set绘制的空间featureplo图改成按照单个基因绘制，可能会有与基因数量过多，导致图画不出来")
}

if ( !is.null(args$help) )	{ print_usage(para) }

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name) || is.null(args$gmt)) {
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
  dir.create(args$outdir, showWarnings = FALSE , recursive = TRUE)
}

if (!file.exists(args$gmt)) {
  cat("Error: gmt file not exists!\n")
  q()
}

if ( is.null(args$reduction)) {
  reduction <- "umap"
}else if (args$reduction == "umap"){
  reduction <- "umap"
}else if (args$reduction == "tsne"){
  reduction <- "tsne"
} else if (args$reduction == "integrated.cca"){
  reduction <- "integrated.cca"
} else if (args$reduction == "integrated.rpca"){
  reduction <- "integrated.rpca"
} else if (args$reduction == "harmony"){
  reduction <- "harmony"
} else if (args$reduction == "umap.harmony"){
  reduction <- "umap.harmony"
} else if (args$reduction == "tsne.harmony"){
  reduction <- "tsne.harmony"
}else if (args$reduction == "integrated.mnn"){
  reduction <- "integrated.mnn"
}else{
  cat("Error: reduction method not exists!\n")
  q()
}

#draw_by_sample <- FALSE
#if ( is.null(args$drawbysample)) {
  #if (args$drawbysample == "T"){
    #draw_by_sample <- TRUE
  #}else if (args$drawbysample == "F"){
    #draw_by_sample <- FALSE
  #}else{
    #cat("Error: drawbysample not exists!\n")
    #q()
  #}
#}



library(Seurat)
library(qusage)
library(ggplot2)
library(scCustomize)
library(patchwork)

gene_plot <- function( pbmc, outdir, gene_set){
    prefix <- paste0(outdir, "/by_gene")
    dir.create(prefix, showWarnings = FALSE)
    for (i in 1:length(gene_set)) {
        name =  names(gene_set[i])
        print( paste0("开始绘制", name, "的图"))
        for (j in 1:length(gene_set[[i]])) {
            if ( ! gene_set[[i]][j] %in% rownames(pbmc) ){
                cat("\n##Error: ", gene_set[[i]][j], "not in rds!\n")
            }else{
                print( paste0("开始绘制", name, "的", gene_set[[i]][j], "的图"))
                a_gene <- gene_set[[i]][j]
                height <- ceiling( sample_number / num_columns) * 4 
                sample_num = length(unique( pbmc$orig.ident))
                print(paste0( "height: ", height))
                tryCatch({
                    pdf(paste0(prefix, "/", name, "_" , a_gene, ".pdf") , width = 16, height = 18)
                    p1 <- DimPlot_scCustom(pbmc, reduction = reduction)
                    p2 <- FeaturePlot_scCustom(pbmc, features= a_gene , reduction = reduction) 
                    p3 <- VlnPlot_scCustom(pbmc, features = a_gene , raster = FALSE)
                    p <- (p1 | p2) / p3 
                    print(p)
                    dev.off()

                    if( sample_num >1 ) {
                        pdf(paste0(prefix, "/", name, "_" , a_gene, "_by_sample.pdf") , width = 9, height = height)
                        p1 <- FeaturePlot_scCustom(pbmc, features= a_gene , reduction = reduction, split.by = "orig.ident" , num_columns = num_columns, repel = TRUE)
                        p2 <- VlnPlot_scCustom(pbmc, features = a_gene , split.by = "orig.ident", raster = FALSE)

                        print(p1,p2)
                        dev.off()
                    }
            }, error = function(e) {
                print(e)
                cat("Error: ", gene_set[[i]][j], "not in rds!\n")
            })
            }
        }

        }
} 


### 绘制所有基因的DotPlot图
gene_dotplot <- function(pbmc, outdir, gene_set, prefix){
    print("开始绘制所有基因的DotPlot图")
    all_gene <- unique(Reduce(union, gene_set))
    pdf(paste0(outdir, "/", prefix, "_", "all_gene_dotplot.pdf") , width = 10, height = 10)
    p1 <- DotPlot_scCustom(pbmc, features = all_gene ,  flip_axes = T,   remove_axis_titles = FALSE) 
    print(p1)
    if(length(gene_set) > 1) {
        prefix <- paste0(outdir, "/by_gene_set") 
        dir.create(prefix, showWarnings = FALSE)
        for (i in 1:length(gene_set)) {
            name =  names(gene_set[i])
            print( paste0("开始绘制", name, "的DotPlot"))
            all_genes <- unique(Reduce(union, gene_set[i]))
            pdf(paste(prefix, "/", name,"_dotplot.pdf",sep=""))
            p1 <- DotPlot_scCustom(pbmc, features = all_genes ,  flip_axes = T,   remove_axis_titles = FALSE)
            print(p1)
            dev.off()
        }
    }
}

gene_featureplot <- function( pbmc, outdir, gene_set, prefix, all_gene ){
    print("开始绘制所有基因的FeaturePlot图")
    pdf(paste0(outdir, "/", prefix, "_", "all_gene_featureplot", ".pdf") , width = 10, height = 10)
    p1 <- DimPlot_scCustom(pbmc, reduction = reduction)
    print(p1)

    gene_number <- length(all_gene)
    num_columns <- min(3, gene_number)
    p2 <- FeaturePlot_scCustom(pbmc, features = all_gene , reduction = reduction,num_columns = num_columns, repel = TRUE)
    print(p2)
    dev.off()
}

gene_vlnplot  <- function( pbmc, outdir, gene_set, prefix , all_gene ){
    print("开始绘制所有基因的VlnPlot图") 
    pdf(paste0( outdir , "/",prefix , "_", "all_gene_vlnplot", ".pdf") , width = 10, height = 10)
    p1 <- Stacked_VlnPlot(pbmc, features = all_gene , x_lab_rotate = TRUE)
    print(p1)
    dev.off()
}


sample_umap_plot <- function( pbmc, outdir  ){
    all_samples <- unique(pbmc$orig.ident)
    print("开始绘制所有样本的FeaturePlot图")
    prefix <- paste0(outdir, "/by_sample") 
    dir.create(prefix, showWarnings = FALSE)

    for (i in all_samples) {
        print(paste0("开始绘制", i, "的FeaturePlot图"))
        pbmc_subset <- subset(pbmc, subset = orig.ident== i)
        prefix_of_sample <- paste0(prefix, "/", i)
        dir.create(prefix_of_sample, showWarnings = FALSE)

        print(paste0("开始绘制", i, "的UMAP图"))
        pdf(paste0(prefix_of_sample, "/", i, "_umap.pdf") , width = 8, height = 6)
        p1 <- DimPlot_scCustom(pbmc_subset, reduction = reduction) + ggtitle( paste0("UMAP of ",  i))
        print(p1)
        dev.off()
        
        print(paste0("开始绘制", i, "的FeaturePlot图"))
        height <- ceiling(length(all_gene) / 3) * 4 
        
        pdf(paste0(prefix_of_sample, "/", i, "_featureplot", ".pdf") , width = 12, height = height)
        p3 <- FeaturePlot_scCustom(pbmc_subset, features = all_gene , reduction = reduction , num_columns = num_columns, repel = TRUE) + ggtitle( paste0("FeaturePlot of ",  i))
        #p <- wrap_plots(p1, p3, ncol = 1)
        print(p3)
        dev.off()
    }
}


gene_heatmap <- function( pbmc, outdir,assay ){
    print("开始绘制每个基因集的HeatMap")
    prefix <- paste0(outdir, "/by_gene_set")
    dir.create(prefix, showWarnings = FALSE)
    for (i in 1:length(gene_set)) {
        name =  names(gene_set[i])
        print( paste0("开始绘制", name, "的HeatMap"))
        pdf(paste0(prefix, "/", name, "_heatmap.pdf") , width = 10, height = 10)
        p1 <- DoHeatmap(pbmc, features = gene_set[[i]],  assay = assay, label = T, group.bar = T, raster = F)
        print(p1)
        dev.off()
    }
}

spatial_plot <- function( pbmc, outdir,gene_set, out_name  ){
	print("开始绘制空间图")
	prefix <- paste0(outdir, "/spatial")
	dir.create(prefix, showWarnings = FALSE)
	print("开始绘制空间的dimplot和空间dimplot图")
	pdf(paste0(prefix, "/", out_name, "_spatial.pdf") , width = 10, height = 10)
	p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, raster = FALSE)
	print(p1)
	p2 <- SpatialDimPlot(pbmc, label = TRUE, label.size = 1, stroke = NA)
	print(p2)
	dev.off()
	print("绘制每个基因空间feaureplot图")
    all_gene <- unique(Reduce(union, gene_set))
    print(paste0("共有基因数量为:",length(all_gene)))
	for (i in 1:length(all_gene)) {
        gene <- all_gene[i]
        print(paste( i, gene, sep=","))
        if(gene %in% rownames(pbmc)){
            pdf(paste0(prefix, "/", out_name, "_", gene, "_spatial.featureplot.pdf") , width = 10, height = 10)
            p3 <- SpatialFeaturePlot(pbmc, features=gene)
            print(p3)
            dev.off()
        }else{
            print(paste0(gene, " not in rds"))
        }
	}
}
########################### 执行部分 ####################################
#读取gene.gmt文件
dir.create(args$outdir, showWarnings = FALSE)

gene_set <- qusage::read.gmt(args$gmt)
all_gene <- unique(Reduce(union, gene_set))

print("读取RDS文件")
pbmc <- readRDS(args$input)
print("读取RDS文件完成")

sample_number <- length(unique(pbmc$orig.ident))
num_columns <- min(3, sample_number)

#按照基因绘图（dimplot, featureplot, vlnplot)，以及在每个样本中的分布图
gene_plot( pbmc, args$outdir, gene_set)
#绘制所有基因的dotplot图
gene_dotplot( pbmc, args$outdir, gene_set, args$name)
### 绘制所有基因的DimPlot图
gene_featureplot( pbmc, args$outdir, gene_set, args$name, all_gene )
### 绘制所有基因的VlnPlot图
gene_vlnplot( pbmc, args$outdir, gene_set, args$name, all_gene )



### 绘制每个样本的UMAP图,如果样本数量超过1绘制
all_samples <- unique(pbmc$orig.ident)
if(length(all_samples) >1 && draw_by_sample ){
    sample_umap_plot( pbmc, args$outdir)
}


assay <- "RNA"
if (!is.null(args$assay)){
  assay <- args$assay
}

### 绘制每个基因集的HeatMap
if(assay != "Spatial"){
    gene_heatmap(pbmc, args$outdir, assay)
}


###绘制空间的图
if(assay == "Spatial"){
    spatial_plot(pbmc, args$outdir, gene_set, args$name)
}




