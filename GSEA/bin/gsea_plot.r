
library('getopt')
para<- matrix(c(
    'help', 'h',    0,  'logical',
    'infile',   'i',    1,  'character',
    'outdir',   'o',    1,  'character',
    'prefix',   'p',    1,  'character'
),byrow=TRUE, ncol=4)

opt <- getopt(para, debug=FALSE)
print_usage <- function(para=NULL) {
    cat(getopt(para, ueasge=TRUE))
    cat("参数说明：
    --infile,-i,[required]输入文件
    --outdir,-o,[required]输出目录
    --prefix,-p,[required]输出文件前缀
    示例：
    Rscript gsea_plot.r -i infile -o outdir -p test
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


## 
data <- read.table( infile, header=T, sep='\t')
for (i in (1:nrow(data))){
    tmp = data[i,]
    id = tmp$ID
    score = tmp$enrichmentScore
    title = paste( id,"_",score,sep="")
    pdf_file = paste(outdir,"/", prefix, title,".gsea.pdf"sep="")
    pdf(pdf_file)
    p1<-gseaplot2(tmp, geneSetID =id,pvalue_table=T,title= title)
    print(p1)
}

