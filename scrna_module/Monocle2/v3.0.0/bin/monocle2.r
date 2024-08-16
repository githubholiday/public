
library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名',
  'ident', 't', 2, 'character', '聚类列名, 默认是seurat_cluster, 选为任意列名' ,
  'root' , 'r', 2, 'character', '根节点，默认是0，可选为ident中某一个名字',
  'draw_column', 'd', 2, 'character', '绘制的列名，默认是seurat_clusters，可选为任意列名,以，分割',
  'assay' , 'a', 2, 'character', '使用的assay，默认是RNA，可选为RNA,ADT,ATAC',
  'only_draw', 'w', 0, 'logical', '是否只绘制聚类图，默认为FALSE',
  'max_components', 'm', 2, 'numeric', '最大的主成分数，默认为2',
  'order_of_ident', 'e', 2, 'character', '图中的顺序, 以，分割'
   ),
  byrow=T,ncol=5
)

args=getopt(command)
if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name)) {
  cat("Usage: Rscript monocle2.r -i input.rds -o outdir -n name\n")
  cat(paste(getopt(command, usage = T), "\n"))
    q() 
}

if (!file.exists(args$input)) {
  cat("Error: input file not exists!\n")
  q()
}

if (!dir.exists(args$outdir)) {
  cat("Warn: outdir not exists!\n")
  dir.create(args$outdir , recursive = TRUE)
}

prefix <- paste0(args$outdir, "/", args$name)


if ( is.null(args$ident)) {
  ident_name <- "seurat_clusters"
}else{
    ident_name <- args$ident
}
print(paste("使用", ident_name))

if ( is.null(args$root)) {
  root <- 0
}else{
    root <- args$root
}
print(paste("使用", root))


if ( is.null(args$draw_column)) {
  draw_column <- c("seurat_clusters")
}else{
    draw_column <- strsplit(args$draw_column, ",")[[1]]
}


if ( is.null(args$assay)) {
  assay <- "RNA"
}else{
    assay <- args$assay
}

if ( is.null(args$only_draw)) {
  only_draw <- FALSE
}else{
    only_draw <- args$only_draw
}

if ( is.null(args$max_components)) {
  max_components <- 2
}else{
    max_components <- args$max_components
}

if ( is.null(args$order_of_ident)) {
  order_of_ident <- c()
}else{
    order_of_ident <- strsplit(args$order_of_ident, ",")[[1]]
}



require(monocle)
require(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
#library(configr)
library(cowplot)

source("/work/share/acuhtwkcu9/taoxiao/05_tool/Stable/Monocle2/v3.0.0/bin/lib/plot_pseudotime_heatmap_get_gene.r")
print("读取RDS"  )



#infile <- "/work/share/acuhtwkcu9/taoxiao/04_Project/P2024042410ZFSUFQ_horse/03_monocle/monocle2/lymphocyte/lymphocyte_0.6.umap.rds"
#pbmc<-readRDS(infile)
#tmp <- immune.combined
#tmp$Barcode <- rownames(tmp@meta.data)
#cluster <- unique(tmp@meta.data$SubGroup)
#meta<-tmp@meta.data
#celllist <- c()
#print(cluster)
#tmp_seurat <- tmp


importCDS2seurat3 <- function (otherCDS , assay="RNA"){
    data <- GetAssayData( otherCDS, assay=assay,slot='counts')
    pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
    fData <- data.frame(gene_short_name = row.names(data),row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)

    expressionFamily <- negbinomial.size()
    expr <- "negbinomial.size"
    monocle_cds <- newCellDataSet(as(data, "sparseMatrix"), phenoData = pd, featureData = fd,expressionFamily = expressionFamily)
    return(monocle_cds)
}

if (only_draw == FALSE){

  pbmc <- readRDS(args$input)
  if ( !root %in% unlist(pbmc[[ ident_name ]] )){
    cat("Error: 根节点不存在，请检查\n")
    q()
  }
  print("转换为monocle2格式")
  HSMM<-importCDS2seurat3(pbmc , assay)
  ## 重新前处理
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  total_cell_number <- ncol(HSMM)
  expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 0.1*total_cell_number))
  #expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
  
  ## 检查ident_name 是否存在
  if (ident_name %in% colnames(pData(HSMM))){
    print(paste0("使用", ident_name, "进行差异分析"))
  }else{
    print(paste0("使用", ident_name, "不存在，请检查"))
    q()
  }

  print(paste0("使用",ident_name, "进行差异分析"))
  diff_test_res <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr = paste0("~", ident_name))
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  HSMM <- setOrderingFilter(HSMM, ordering_genes)
  HSMM <- reduceDimension(HSMM, max_components = max_components,  method = 'DDRTree')
  HSMM <- orderCells(HSMM)

  GM_state <- function(cds , ident , root  ){
    if (length(unique(pData(cds)$State)) > 1){
      T0_counts <- table(pData(cds)$State, pData(cds)[[ident_name ]] )[, root ]
      return(as.numeric(names(T0_counts)[which
            (T0_counts == max(T0_counts))]))
    } else {
      return (1)
    }
  }

  HSMM <- orderCells(HSMM, root_state = GM_state(HSMM , ident_name , root ))
  saveRDS(HSMM, paste0(prefix, ".rds"))

}else{
  HSMM <- readRDS(args$input)
  #expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
	total_cell_number <- ncol(HSMM)
	expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 0.1*total_cell_number))

  print("只绘制图")
}


print("绘制轨迹DDTree")
pdf(paste0(prefix, ".trajectory.pdf"))
p1 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime" )
if (length(order_of_ident) > 0){
  pData(HSMM)[[ident_name ]] <- factor(pData(HSMM)[[ident_name ]] ,  levels = order_of_ident)
  p2 <- plot_cell_trajectory(HSMM, color_by = ident_name) + scale_color_discrete(limits = order_of_ident)+scale_colour_manual( values=c("#f2a7da","#71b7ed","#F9B562"))
}else{
  p2 <- plot_cell_trajectory(HSMM, color_by = ident_name)

}
p3 <- plot_cell_trajectory(HSMM, color_by = "State")
print(p1)
print(p2)
print(p3)
for ( i in draw_column){
  p4 <- plot_cell_trajectory(HSMM, color_by = i)
  print(p4)
}
dev.off()


print("统计state和ident的交集情况")
table_result <- table(pData(HSMM)$State , pData(HSMM)[[ident_name]])
write.table(table_result , paste0(prefix,".State_ident.xls"),sep="\t",quote=F,row.names=T)




print("对State进行差异分析")
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],  fullModelFormulaStr = "~State")
diff_test_res_short <-  diff_test_res[,c("gene_short_name", "pval", "qval")]
write.table(diff_test_res_short , paste0(prefix,".DEG_State_genes.xls"),sep="\t",quote=F,row.names=F)
### 选取top10的DEG，绘制图片
top_genes <- diff_test_res %>%
  arrange(qval) %>%
  head(10) %>%
  pull(gene_short_name)
hsmm_subset <- HSMM[top_genes,]
pdf(paste0(prefix,".DEG_State_genes.pdf"))
plot_genes_jitter(hsmm_subset,
                  grouping="State",
                  color_by="State",
                  nrow=NULL,
                  ncol=3,
                  plot_trend=TRUE)
dev.off()


print("寻找随伪时间变化的基因")
diff_test_res <- differentialGeneTest(HSMM[expressed_genes, ],  fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res_short <-  diff_test_res[,c("gene_short_name", "pval", "qval")]
write.table(diff_test_res_short , paste0("DEG_Pseudotime_genes.xls"),sep="\t",quote=F,row.names=F)

top_genes <- diff_test_res %>%
  arrange(qval) %>%
  head(10) %>%
  pull(gene_short_name)

hsmm_subset <- HSMM[top_genes,]
pdf(paste0(prefix,".DEG_Pseudotime_genes.pdf"))
p1 <- plot_genes_in_pseudotime(hsmm_subset, color_by = "Pseudotime")
print(p1)
p2 <- plot_genes_in_pseudotime(hsmm_subset, color_by = "orig.ident")
print(p2)
p3 <- plot_genes_in_pseudotime(hsmm_subset, color_by = "State")
print(p3)
dev.off()



sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
pdf(paste0(prefix,".DEG_Pseudotime_genes_heatmap.pdf"))
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                num_clusters = 3,
                cores = 1,
                show_rownames = T)
dev.off()

annotation_row <- plot_pseudotime_heatmap2(HSMM[sig_gene_names,],num_clusters = 3,cores = 1,show_rownames = T)
gene_outfile = paste0(prefix, ".heatmap.gene.cluster.xls")
write.table( annotation_row, gene_outfile, sep="\t", quote=FALSE, row.names=FALSE)

## 分析单细胞轨迹中的分支
print("分析单细胞轨迹中的分支")
tryCatch({
  BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 1)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  write.table(BEAM_res , paste0(prefix,".BEAM_genes.xls"),sep="\t",quote=F,row.names=F)
  pdf(paste0(prefix,".BEAM_genes_heatmap.pdf"))
  plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,
                                            qval < 1e-4)),],
                                            branch_point = 1,
                                            num_clusters = 4,
                                            cores = 1,
                                            use_gene_short_name = T,
                                            show_rownames = T)
  dev.off()

  top_genes <- BEAM_res %>%
    arrange(qval) %>%
    head(10) %>%
    pull(gene_short_name)


  #lung_genes <- row.names(subset(fData(HSMM),  gene_short_name %in%  top_genes  ))
  pdf(paste0(prefix,".BEAM_genes_branch.pdf"))
  plot_genes_branched_pseudotime(HSMM[top_genes,],
                        branch_point = 1,
                        color_by = "Time",
                        ncol = 1)
  dev.off()
},
error = function(e) {
  print(e)
} 
)
