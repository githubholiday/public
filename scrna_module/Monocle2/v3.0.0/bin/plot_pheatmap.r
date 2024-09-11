#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
#setwd("F:\\plot")
#最初用途是为了调整monocle2中heatmap图的分组顺序，也可以用于所有的pheatmap的分组顺序调整
#使用跑完monocle2的rds中输出heatmap图中的matrix文件，手动调整重新聚类后的cluster在热图中的顺序
#顺序处需要手动调整，暂时无法通过参数实现

library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '跑完monocle2的rds文件',
  'outfile', 'o', 1, 'character', '输出的pdf文件全路径'
   ),
  byrow=T,ncol=5
)

args=getopt(command)
if (!is.null(args$help) || is.null(args$input) || is.null(args$outfile)) {
  cat("Usage: Rscript monocle2.r -i input.rds -o outfile \n")
  cat(paste(getopt(command, usage = T), "\n"))
    q() 
}


require(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
library(cowplot)
library(pheatmap)


#颜色操作
table.ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1)
{
  x <- seq(0, 1, length.out = n)
  y <- rep(0, length(x))
  sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
  sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
  y[sill.min:sill.max] <- 1
  base.min <- round((n - 1) * (mid - base / 2)) + 1
  base.max <- round((n - 1) * (mid + base / 2)) + 1
  xi <- base.min:sill.min
  yi <- seq(0, 1, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  xi <- sill.max:base.max
  yi <- seq(1, 0, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  height * y
}

rgb.tables <- function(n,
        red = c(0.75, 0.25, 1),
        green = c(0.5, 0.25, 1),
        blue = c(0.25, 0.25, 1))
{
  rr <- do.call("table.ramp", as.list(c(n, red)))
  gr <- do.call("table.ramp", as.list(c(n, green)))
  br <- do.call("table.ramp", as.list(c(n, blue)))
  rgb(rr, gr, br)
}

matlab.like2 <- function(n)
  rgb.tables(n,
   red = c(0.8, 0.2, 1),
   green = c(0.5, 0.4, 0.8),
   blue = c(0.2, 0.2, 1))

blue2green2red <- matlab.like2
## 颜色定义
bks <- seq(-3.1,3.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)
### 颜色操作


#####################################################################
#定义输入和输出文件
infile <- args$input
outpdf <- args$outfile
rawpdf <- paste0(outpdf,".raw.pdf") #调整顺序前的图
#读取输入文件
m <- read.table(infile,header=TRUE)
#处理输入数据
heatmap_matrix <- m
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1

#定义参数
cluster_rows=TRUE
hclust_method = "ward.D2" #聚类方法
num_clusters = 3 #最终聚类数量,可以通过修改，生成不同cluster数量
use_gene_short_name = TRUE
ph <- pheatmap(heatmap_matrix, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=cluster_rows, 
               show_rownames=F, 
               show_colnames=F, 
               clustering_distance_rows=row_dist,
               clustering_method = hclust_method,
               cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols)

#根据上面聚类的结果获取基因的聚类情况
annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))

colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))

ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                   useRaster = T,
                   cluster_cols = FALSE, 
                   cluster_rows = cluster_rows, 
                   show_rownames=T, 
                   show_colnames=F, 
                   #scale="row",
                   clustering_distance_rows=row_dist, #row_dist
                   clustering_method = hclust_method, #ward.D2
                   cutree_rows=num_clusters,
                   # cutree_cols = 2,
                   annotation_row=annotation_row,
                   annotation_col=NA,
                   treeheight_row = 20, 
                   breaks=bks,
                   fontsize = 6,
                   color=hmcols, 
                   border_color = NA,
                   silent=TRUE,
                   filename=rawpdf
)

######################### 从这里开始要手动调图 #######################
#基于原图的顺序，调整分组的顺序
old_order <- ph_res$tree_row$order #现在的顺序应该是 2 1 3 ,要调整为312
#看一下原图是什么顺序 计算一下，调整顺序，输出到manul_ordel里
#查看总的数量，以及各个cluster的数量
dim(annotation_row) #5691,从1开始计数
length(annotation_row[annotation_row[,1]=='1',]) #1952 3740:5691
length(annotation_row[annotation_row[,1]=='2',]) #1421 2319:3739
length(annotation_row[annotation_row[,1]=='3',]) #2318 1:2318
2318+1421
#rownames(heatmap_matrix)[old_order]
manul_order <- old_order[c(1:2318,3740:5691,2319:3739)]
#rownames(heatmap_matrix)[manul_order]

#对数据进行聚类，并且重新排序
hclust1 <- hclust(row_dist,method=hclust_method)
#hclust1[1:5]
denf = reorder(as.dendrogram(hclust1), wts=order(manul_order, rownames(heatmap_matrix)))
row_cluster <- as.hclust(denf,method=hclust_method)

ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                   useRaster = T,
                   cluster_cols = FALSE, 
                   cluster_rows = row_cluster, 
                   show_rownames=T, 
                   show_colnames=F, 
                   #scale="row",
                   clustering_distance_rows=row_dist, #row_dist
                   #clustering_method = hclust_method, #ward.D2
                   cutree_rows=num_clusters,
                   # cutree_cols = 2,
                   annotation_row=annotation_row,
                   annotation_col=NA,
                   treeheight_row = 20, 
                   breaks=bks,
                   fontsize = 4,
                   color=hmcols, 
                   border_color = NA,
                   silent=TRUE,
                   filename=outpdf
)


