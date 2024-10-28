library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'prefix' , 'p', 1, 'character', '输出的文件名',
  'name' , 'n', 1, 'character', '使用的细胞群对象名，一般为seurat_clusters'
),
  byrow=T,ncol=5
)


## 读取参数
args=getopt(command)

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

library(Seurat)
library(qusage)
library(ggplot2)
library(scCustomize)
library(patchwork)

print("读取rds文件")
marrow <- readRDS( args$input)

name <- args$name
table(marrow@meta.data$name)

#获取基因
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#细胞周期计算
marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# 不同name中各细胞周期的数量,并输出
marrow@meta.data$Phase <- factor(marrow@meta.data$Phase,levels=c("G1","S","G2M"))
table(marrow@meta.data$Phase)
tt <- table(marrow@meta.data$Phase,marrow[[name]])
tt <- cbind(Phase=rownames(tt),tt)
stat_file <- paste0(args$outdir, "/", args$prefix, ".cellcycle.xls")
write.table(tt,stat_file,sep="\t",quote=FALSE, row.names=FALSE)

#绘制得分的featureplot
feature_pdf <- paste(args$outdir, "/", args$prefix, ".cellcycle.featureplot.pdf")
pdf(feature_pdf)
FeaturePlot(marrow,features=c("S.Score"))
FeaturePlot(marrow,features=c("G2M.Score"))
dev.off()

#绘制不同cluster中各细胞周期的比例
## 计算细胞比例
Cellratio <- prop.table(table(marrow[["Phase"]],marrow[[name]]))
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("CellCycle",name,"Freq")
head(Cellratio)

## 绘制堆叠图
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

stack_pdf <- paste(args$outdir, "/", args$prefix, ".cellcycle.stackplot.pdf")
pdf(stack_pdf)
ggplot(Cellratio) + 
#geom_bar(stat = "identity",position="fill")+
    geom_bar(aes(x =name, y= Freq, fill = CellCycle),stat = "identity",position="fill",width = 0.7,size = 0.5,colour = '#222222')+ 
    theme_classic() +
    labs(x='Cluster',y = 'Ratio')+
  #scale_fill_manual(values = allcolour)+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()