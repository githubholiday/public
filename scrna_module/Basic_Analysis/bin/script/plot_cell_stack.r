library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'prefix' , 'n', 1, 'character', '输出的文件名'
),
  byrow=T,ncol=5
)


## 读取参数
args=getopt(command)

print_usage <- function(para=NULL){
    cat("参数说明：
    -i:rds文件
    -o:输出目录，可以不存在
    -p:输出前缀

    使用实例：
    /public/software/apps/singularity/3.7.3/bin/singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/liutao/sif/sif/scRNA/seurat5/seurat5_sccustomize.sif Rscript plot_cell_stack.r -i rds -o outdir -p sample ")
}

if ( !is.null(args$help) )	{ print_usage(para) }

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$prefix)) {
  cat("Usage: Rscript cluster_umap -i input.rds -o outdir -p prefix\n")
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


library(ggplot2)
library(Seurat)


#判断两两指标的细胞数量统计
#计算细胞比例
plot_stack <- func( marrow, value1, value2, name="CellType" ){
    Cellratio <- prop.table(table(marrow[[value1]],marrow[[value2]]))
    Cellratio <- as.data.frame(Cellratio)
    #默认是Var1,Var2,Freq，进行重命名
    colnames(Cellratio) <- c(value1,value2,"Freq")

    prefix <- paste(args$outdir, "/", args$prefix, sep="")

    pdf(paste( prefix, "_", name, "_stackplot.pdf",sep=""))

    ggplot(Cellratio) + 
    geom_bar(aes(x = value1, y= Freq, fill = value2),stat = "identity",position="fill",width = 0.7,size = 0.5,colour = '#222222')+ 
    theme_classic() +
    labs(x='Cluster',y = 'Ratio')+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
    dev.off()

}


marrow <- readRDS( args$input )
#绘制样本的图
plot_stack(marrow, "orig.ident", 'seurat_clusters')
plot_stack(marrow, "Group", 'seurat_clusters')
plot_stack(marrow, "orig.ident", 'CellType')
plot_stack(marrow, "Group", 'CellType')

