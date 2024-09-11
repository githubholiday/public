#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
#setwd("F:\\plot")
#最初用途是为了调整monocle2中heatmap图的分组顺序，也可以用于所有的pheatmap的分组顺序调整

#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
#setwd("F:\\plot")
#最初用途是为了调整monocle2中heatmap图的分组顺序，也可以用于所有的pheatmap的分组顺序调整
#使用monocle2跑完rds中导出的matrix文件，使用pheatmap进行聚类（默认聚成3类）
#将重新聚类后的每cluster的基因输出


library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', 'monocle2中导出的matrix文件',
  'outpre', 'o', 1, 'character', '输出的前缀'
   ),
  byrow=T,ncol=5
)

args=getopt(command)
if (!is.null(args$help) || is.null(args$input) || is.null(args$outpre)) {
  cat("Usage: Rscript monocle2.r -i input.rds -o outdir/name \n")
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
#读取输入文件
m <- read.table(args$infile,header=TRUE)
#处理输入数据
heatmap_matrix <- m
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1

#定义参数
cluster_rows=TRUE
hclust_method = "ward.D2" #聚类方法
num_clusters = 3 #最终聚类数量
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

cluster <-  unique(annotation_row[,"Cluster"])
annotation_row <- cbind("Gene"=rownames(annotation_row),annotation_row)
for (i in cluster){
  outfile = paste0(args$outpre,i,".xls")
  gene_name = annotation_row[annotation_row["Cluster"]==i,]
  write.table(gene_name, outfile, quote=FALSE, sep="\t",row.names=FALSE)
  }
