library('getopt')
para<- matrix(c(
        'help', 'h',    0,      "logical", "帮助文档",
        'infile',       'i',    1,      "character","输出的rds文件",
        'outdir',       'o',    1,      "character","结果输出目录，如果不存在会创建",
        'name',         'n',    1,      "character","输出结果文件的前缀",
        'Ident',        'I',    1,      "character","分析使用的标签，如 seurat_clusters",
        'assay',         'a',    1,      "character","数据类型，是转录组数据还是空间数据，转录组数据使用RNA，空间数据使用Spatial"
),byrow=TRUE,ncol=5)
args <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
    cat(getopt(para,usage=TRUE))
    cat("
    功能：使用cellchat对单细胞转录组数据进行细胞通讯分析
    
    参数说明：
    --infile, -i :[required] 输入的rds文件
    --outdir, -o :[required] 输出目录，如果不存在，则会创建
    --name,   -n :[required] 输出文件前缀信息 
    --Ident,  -I :[required] 分析使用的标签信息，如 seurat_clusters 
    --assays, -a :[alternative] 数据类型，可选 [RNA, Spatial], RNA=单细胞转录组数据，Spatial=空间转录组数据，默认RNA 

    示例：
    singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9 /work/share/acuhtwkcu9/taoxiao/07_sif/scRNA/cellchat.sif Rscript  CellChat_seurat5.r -i rds -o oudir -n sample -I CellType\n")
    q(status=1)
}

if (!is.null(args$help) || is.null(args$infile) || is.null(args$outdir) || is.null(args$name) || is.null(args$Ident)) {
    print_usage(para) 
}

if(is.null(args$assay )){
    data_assays <- args$assay
}else{
    data_assays <- "RNA"
}
#定义输入参数
rds <- args$infile
outdir <- args$outdir
prefix <- args$name
out_pre = paste( outdir, prefix , sep="/" )
Ident <- args$Ident

#创建输出目录
if(dir.exists( outdir )){
    print(paste(outdir," already exists!"))
}else{
    dir.create( outdir )
}

#导入需要的包
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
## 读取数据
pbmc<-readRDS(rds)
#pbmc<- subset (pbmc , subset = group == "POEMS")
print("读取数据完成")

## 创建CellChat对象
cellchat <- createCellChat(object =  pbmc, group.by = Ident, assay = data_assays )

## 导入数据库--人数据库
#CellChatDB <- CellChatDB.human
#cellchat@DB <- CellChatDB
## 导入数据库--小鼠
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB

## 计算细胞间通讯
#future::plan("multisession", workers = 4) # do parallel
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


cellchat <- computeCommunProb(cellchat, type = "triMean")


cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

## 导出结果

df.net <- subsetCommunication(cellchat)

output_data_frame <- function(data , output){
    data<-cbind("ID" = rownames(data),data)
    write.table(data,file = output,sep = "\t",quote = FALSE,row.names = FALSE)
}

output_data_frame(df.net, paste(out_pre,"cellchat_weight.txt",sep="."))


groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
pdf(paste(out_pre,"cellchat_all_interaction_count.pdf",sep="."))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf(paste(out_pre,"cellchat_all_interaction_weight.pdf",sep="."))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
#par(mfrow = c(3,4))
pdf(paste(out_pre,"cellchat_weight.pdf",sep="."))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  print(rownames(mat)[i])
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = paste0( "interaction weights of " , rownames(mat)[i]))
}
dev.off()


mat <- cellchat@net$count
#par(mfrow = c(3,4))
pdf(paste(out_pre,"cellchat_count.pdf",sep="."))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  print(rownames(mat)[i])
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = paste0( "interaction count of " , rownames(mat)[i]   ))
}
dev.off()


pdf(paste(out_pre,"cellchat_weight_heatmap.pdf",sep="."))
netVisual_heatmap( cellchat , measure = "weight")
dev.off()

pdf(paste(out_pre,"cellchat_count_heatmap.pdf",sep="."))
netVisual_heatmap( cellchat , measure = "count")
dev.off()

saveRDS(cellchat, file = paste(out_pre,"cellchat.rds",sep="."))

