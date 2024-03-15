
library('getopt')
para<- matrix(c(
    'help', 'h',    0,  'logical',
    'rds',  'r',      1, 'character',
    'idfile',   'i',    1,  'character',
    'outdir',   'o',    1,  'character',
    'prefix',   'p',    1,  'character',
    'header',   'header', 1,  'logical'
),byrow=TRUE, ncol=4)

opt <- getopt(para, debug=FALSE)
print_usage <- function(para=NULL) {
    cat(getopt(para, ueasge=TRUE))
    cat("
【功能说明】
    对GSEA的结果进行经典图绘制,输入为一列ID即可。图中展示NES,pvalue, p.adjust三列信息
【参数说明】
    --rds,-r,[required]输入文件,gsea结果的rds文件
    --idfile,-i,[required]需要绘图的ID列表输入文件（可以包含表头，也可以不包含）
    --outdir,-o,[required]输出目录
    --prefix,-p,[required]输出文件前缀
    --header,-header,[optional]是否包含标题,默认不包含
    --help,-h,[optional]帮助文档
【使用示例】
    Rscript gsea_plot.r -r gsea.rds -i infile -o outdir -p test -header False
    ")
    q(status=1)
}

#对输入进行检查
if ( !is.null(opt$help) )>--{ print_usage(para) }
if ( is.null(opt$infile) )>--{ cat("Please input the data file\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )>--{ cat("Please input the outdir of result\n\n") ; print_usage(para)}
if ( is.null(opt$prefix) )>--{ cat("Please input the prefix of result\n\n") ; print_usage(para)}

#对输出进行赋值处理
infile <- para$infile
outdir <- para$outdir
prefix <- para$prefix

library(enrichplot)
library(ggplot2)
library(ggpp)


## 
rds_data <- readRDS(rds_data)
ids_data <- read.table( infile, header=T, sep='\t')
outpre <- paste( outdir, prefix, sep='/')
GSEA_plot( rds_data, outpre, ids )

GSEA_plot<-function(gsea_res,out_pre, ids ){
    library(ggpp)
    #如果想将所有id的画在一个图中，不循环即可
    for (k in (1:length(ids))){
        geneSetID = ids[k]
        filename<-paste(out_pre, geneSetID,'GSEA_score.pdf',sep='_')
        pdf(filename, w=12, h=8)
        
        pd <- gsea_res[geneSetID, c( "NES","pvalue", "p.adjust")]
        rownames(pd) <- NULL
        for (i in seq_len(ncol(pd))){pd[, i] <- format(pd[, i], digits = 4)}
        pd_t <- data.frame(t(pd))
        data<-cbind("Type" = rownames(pd_t),pd_t)
        colnames(data) <- c("Type","Value")
        pd <- data
        
        #绘图
        p <- gseaplot2(gsea_res, geneSetID = geneSetID)

        p[[1]] <- p[[1]]+
            annotate("table", x = 12000, y = 0.4, label = pd,
            size=3,table.theme = ttheme_gtminimal
            )+
            theme(plot.title = element_text(size = 5),
                legend.position = "top",
                legend.direction = "vertical"
            )
        print(p)
        dev.off()
        }
}

