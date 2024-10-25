library('getopt')

para<- matrix(c(
 'help', 'h', 0, 'logical',
 'rds',  'i', 1, 'character',
 'outdir',  'o',  1, 'character',
 'prefix',  'p',  1, 'character',
 'Ident',  'I',  1, 'character'
),byrow=TRUE, ncol=4)
  
opt <- getopt(para, debug=FALSE)
print_usage <- function(para=NULL) {
    cat(getopt(para, ueasge=TRUE))
    cat("
【功能说明】
    对单个样本做差异基因分析，使用FindAllMarkers函数对cluster做差异基因分析,并输出上下调和是否显著
【参数说明】
    --rds,-r,[required]输入的rds文件
    --outdir,-o,[required]输出目录
    --prefix,-p,[required]输出文件前缀
    --help,-h,[optional]帮助文档
【使用示例】
    Rscript gsea_plot.r -r infile -o outdir -p test")
     q(status=1)
}

#对输入进行检查
if ( !is.null(opt$help) )   { print_usage(para) }
if ( is.null(opt$infile) )  { cat("Please input the data file\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )  { cat("Please input the outdir of result\n\n") ; print_usage(para)}
if ( is.null(opt$prefix) )  { cat("Please input the prefix of result\n\n") ; print_usage(para)}

  
#对输入进行赋值处理
infile <- para$infile
outdir <- para$outdir
prefix <- para$prefix
 
library(Seurat)
library(dplyr)

rds <- readRDS(opt$rds)
Idents(rds) <- args$Ident
all.markers <- FindAllMarkers( rds,logfc.threshold = 0.1, test.use = "wilcox",min.pct = 0.01 )
#增加基因列
all.markers <- cbind("Gene"=rownames(all.markers),all.markers)

#判断上下调和是否显著
all.markers <- all.markers %>% mutate(up_down = ifelse(avg_log2FC > 0, "Up", "Down"),significant = ifelse(p_val_adj < 0.05, "yes", "no"))
out_pre <- paste(opt$outdir, $prefix, sep="/")

for (i in unique(sto@meta.data$seurat_clusters)){
	markers <- all.markers[ all.markers$cluster==i,]
	outfile <- paste(out_pre, i, "marker.xls",sep="_")
	write.table(markers, outfile, sep="\t", quote=FALSE,row.names = FALSE )
    }
	
