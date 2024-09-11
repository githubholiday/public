#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
#setwd("F:\\plot")
#最初用途是为了调整monocle2中heatmap图的分组顺序，也可以用于所有的pheatmap的分组顺序调整
#从跑完monocle2的rds中输出heatmap图中的matrix文件
#文件第一行为cluster(100个)，第一列为基因，值为基因在不同cluster中的表达值

library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '跑完monocle2的rds文件',
  'outfile', 'o', 1, 'character', '输出的文件'
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
library(monocle)
#定义输入和输出
monocle2_rds <- args$input
#outfile <- paste(args.outdir, "/", "matrix.xls", sep="" )

HSMM <- readRDS(monocle2_rds)

total_cell_number <- ncol(HSMM)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 0.1*total_cell_number))
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],  fullModelFormulaStr = "~State")
diff_test_res_short <-  diff_test_res[,c("gene_short_name", "pval", "qval")]

diff_test_res <- differentialGeneTest(HSMM[expressed_genes, ],  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
cds_subset <- HSMM[sig_gene_names,]


#定义参数：plot_pseudotime_heatmap函数中国的参数
num_clusters = 3 #默认是6，咱们脚本里是3
num_clusters <- min(num_clusters, nrow(cds_subset))
pseudocount <- 1
trend_formula = '~sm.ns(Pseudotime, df=3)'
scale_max=3
scale_min=-3
cores=1

###### 

newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), max(pData(cds_subset)$Pseudotime),length.out = 100)) 

m <- genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,  relative_expr = T, new_data = newdata)

m <- m[!apply(m,1,sum)==0,]
m = log10(m+pseudocount)

m=m[!apply(m,1,sd)==0,]
m=Matrix::t(scale(Matrix::t(m),center=TRUE))
m=m[is.na(row.names(m)) == FALSE,]
m[is.nan(m)] = 0
m[m>scale_max] = scale_max
m[m<scale_min] = scale_min

write.table(m, args$outfile, quote=FALSE,sep="\t")
