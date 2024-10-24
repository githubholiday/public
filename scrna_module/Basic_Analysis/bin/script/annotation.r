library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名',
  'database' , 'd', 1, 'character', '数据库的名称，可选为1: BlueprintEncodeData，2:DatabaseImmuneCellExpressionData , 3:HumanPrimaryCellAtlasData, 4:ImmGenData , 5:MonacoImmuneData , 6:MouseRNAseqData,  7:NovershternHematopoieticData'
),
  byrow=T,ncol=5
)

## 读取参数
args=getopt(command)

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name)) {
  cat("Usage: Rscript cluster_umap -i input.rds -o outdir -n name\n")
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
print("输出目录已创建")


library(Seurat)
library(SingleR)
library(gplots) 
library(knitr)

library(celldex)


pbmc <- readRDS(args$input)


if ( is.null(args$database)) {
  ref <- readRDS("/dabase/R_data/celldex/BlueprintEncode.rdata")
  database_name <- "BlueprintEncodeData"
}else if (args$database == "1"){
  ref <- readRDS("/dabase/R_data/celldex/BlueprintEncode.rdata")
  database_name <- "BlueprintEncodeData"
}else if (args$database == "2"){
  ref <- readRDS("/dabase/R_data/celldex/DatabaseImmuneCellExpression.rdata")
  database_name <- "DatabaseImmuneCellExpressionData"
}else if (args$database == "3"){
  ref <- readRDS("/dabase/R_data/celldex/HumanPrimaryCellAtlas.rdata")
  database_name <- "HumanPrimaryCellAtlasData"
}else if (args$database == "4"){
  ref <- readRDS("/dabase/R_data/celldex/ImmGen.rdata")
  database_name <- "ImmGenData"
}else if (args$database == "5"){
  ref <- readRDS("/dabase/R_data/celldex/MonacoImmune.rdata")
  database_name <- "MonacoImmuneData"
}else if (args$database == "6"){
  ref <- readRDS("/dabase/R_data/celldex/MouseRNAseq.rdata")
  database_name <- "MouseRNAseqData"
}else if (args$database == "7"){
  ref <- readRDS("/dabase/R_data/celldex/NovershternHematopoietic.rdata")
  database_name <- "NovershternHematopoieticData"
}else{
  cat("Error: database not exists!\n")
  q()
}

pbmc.data<-GetAssayData(pbmc,layer="data")

clusters <- pbmc$seurat_clusters # 获取现在的细胞类群名字
table(clusters)

prefix <- paste0(args$outdir, "/", args$name, "_", database_name , "_" )

## 对每个细胞进行注释
pbmc.singler.percell<-SingleR(test = pbmc.data, ref = ref, labels = ref$label.main )

prefix_percell <- paste0(prefix, "percell_" )

pdf(file = paste(prefix_percell, "score_heatmap.pdf", sep = "") , width = 10, height = 10)
p1 <- plotScoreHeatmap(pbmc.singler.percell, clusters = clusters)
print(p1)
p2<- plotScoreHeatmap(pbmc.singler.percell)
print(p2)
dev.off()

singleR_percell_annotations<-table(pbmc.singler.percell$labels , pbmc$seurat_clusters )

show_table_function_result <- function(table_result , prefix){
  library(ggpubr)
  write.csv(as.matrix(table_result) , file=paste0(prefix,".csv"))
  
  pdf(paste0(prefix,".pdf"),w=8,h=6)
  p <- ggballoonplot(as.data.frame(table_result) , fill = "value", color = "lightgray", size = 10, show.label = TRUE) +  gradient_fill(c("blue", "white", "red")) 
  print(p)
  dev.off()
}

show_table_function_result(singleR_percell_annotations , paste0(prefix_percell , "stat"))



write.table(as.matrix(pbmc.singler.percell) , file=paste(prefix_percell, "annotation.xls", sep = "") , sep = "\t")


### 对每个细胞类群进行注释
pbmc.singler.cluster <-SingleR(test = pbmc.data, ref = ref, labels = ref$label.main  , clusters = clusters)

prefix_cluster <- paste0(prefix, "cluster_" )
pdf(file = paste(prefix_cluster, "score_heatmap.pdf", sep = "") , width = 10, height = 10)
p1 <- plotScoreHeatmap(pbmc.singler.cluster, clusters = clusters)
print(p1)
p2<- plotScoreHeatmap(pbmc.singler.cluster)
print(p2)
dev.off()

write.table(as.matrix(pbmc.singler.cluster) , file=paste(prefix_cluster, "annotation.xls", sep = "") , sep = "\t")





#new.clusterID <- pred.immune$labels
#names(new.clusterID) <- levels(immune.combined)
#immune.combined.new <- RenameIdents(immune.combined,new.clusterID)



#plotScoreHeatmap(pbmc.singler)
#plotDeltaDistribution(pbmc.singler, ncol = 3)

#summary(is.na(pbmc.singler$pruned.labels))



#saveRDS(pbmc,file="pbmc.unintegrated.rds")
#genes_to_check = c("MITF","TYR","SOX10","S100B","PTPRC","CD3E","GZMB","GZMA","CD79A","MZB1","SDC1","LYZ",'PECAM1',"COL6A2")
#p<-DotPlot(pbmc,features = genes_to_check,assay="RNA")

