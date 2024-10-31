library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名' ,
  'min.pct', 'p', 2, 'numeric', '最小百分比，默认为0.01',
  'logfc.threshold', 'l', 2, 'numeric', 'logfc阈值，默认为0.1',
  'test.use', 't', 2, 'character', '检验方法，默认为wilcox，可选为wilcox或roc等',
  'Ident' , "I" , 2 , 'character' , '聚类的标签，默认为harmony_clusters， 可选为任意列名'
),
  byrow=T,ncol=5)

## 读取参数
args=getopt(command)

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name) ){
  cat("Usage: Rscript Pseudobulk.r -i input.rds -o outdir -n name\n")
  cat(paste(getopt(command, usage = T), "\n"))
    q() 
}

if (!file.exists(args$input)) {
  cat("Error: input file not exists!\n")
  q()
}

if (!dir.exists(args$outdir)) {
  dir.create(args$outdir , recursive = T)
  q()
}
print(paste0("outdir:",args$outdir))

if ( is.null(args$min.pct)) {
  min.pct <- 0.01
}else{
    min.pct <- args$min.pct
    }

print(paste0("min.pct:",min.pct))


if ( is.null(args$logfc.threshold)) {
    logfc.threshold <- 0.1
}else{
    logfc.threshold <- args$logfc.threshold
    }

print(paste0("logfc.threshold:",logfc.threshold))

if ( is.null(args$test.use)) {
    test.use <- "wilcox"
}else{
    test.use <- args$test.use
    }


print(paste0("test.use:",test.use))
library(dplyr)

if ( is.null(args$Ident)) {
    Ident.col <- "harmony_clusters"
}else{
    Ident.col <- args$Ident
}

print(paste0("Ident:",Ident.col))


output_data_frame <- function(data , output){
    data<-cbind("Gene"= rownames(data),data)
    write.table(data,file = output,sep = "\t",quote = FALSE,row.names = FALSE)
}


library(Seurat)

pbmc <- readRDS(args$input)

print(paste0("检查meta信息中是否有",Ident.col,"信息"))
if ( !Ident.col %in% colnames(pbmc@meta.data)){
  cat("Error: meta信息中没有",Ident.col,"信息!, 请检查\n")
  q()
}else{
  Idents(pbmc) = Ident.col
}


all.markers <- FindAllMarkers(object = pbmc , only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold, test.use = test.use)
#all.markers <- all.markers[,c(7,6,1,2,3,4,5)]
head(all.markers)
#增加基因列
#all.markers <- cbind("Gene"=all.markers[,"gene"],all.markers)
names(all.markers)[names(all.markers)=="gene"] <- "Gene"
all.markers <- all.markers[,c(7,1:6)]

#判断上下调以及是否显著
all.markers <- all.markers %>% mutate(up_down = ifelse(avg_log2FC > 0, "Up", "Down"),significant = ifelse(p_val_adj < 0.05, ifelse(abs(avg_log2FC) > 0.25, "yes", "no"), "no"))
#输出总的all_marker表
write.table(all.markers,file = paste0(args$outdir,"/",args$name,".markers.xls") ,sep = "\t",quote = FALSE,row.names = FALSE)

write.table(pbmc$harmony_clusters,file = paste0(args$outdir,"/",args$name,".cluster.xls") ,sep = "\t",quote = FALSE ,col.names = FALSE)

#按照每个cluster输出marker
all_cluster <- sort(unique( unlist(pbmc[[Ident.col]]) ))
for (i in all_cluster){
	outfile <- paste0(args$outdir,"/",args$name,".cluster",i,".markers.xls")
	cluster_data <- all.markers %>% filter( all.markers[,"cluster"] %in% i)
	write.table(cluster_data,file=outfile,sep = "\t",quote = FALSE,row.names = FALSE)
}
print("完成")