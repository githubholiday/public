library('getopt')
para<- matrix(c(
    "help"        ,"h",     0,    "logical",
    "rds"         ,"i",     2,    "character",
    "db"          ,"d",     2,    "character",
    "species"     ,"s",     2,    "character",
    "group"       ,"g",     1,    "character",
    "ident"       ,"c",     1,    "character",
    "outdir"      ,"o",     2,    "character",
    "multiprocess","m",     1,    "numeric"),
    byrow=TRUE,ncol=4)
#===========================================================
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
    cat(getopt(para,usage=TRUE))
    cat("
功能说明:仅针对人和小鼠的物种做细胞通讯分析,注意选择正确的物种
    human(1939对229种pathwayL-R互作); mouse(2019对229种pathwayL-R互作)
    Secreted Signaling,ECM-Receptor,Cell-Cell Contact
    human:1199/421/319; mouse:1209/432/378
参数说明:
    --help:帮助文档
    --rds:[必需]rds文件(单细胞转录组)
    --db:[必需]default Secreted Signaling],Secreted Signaling,ECM-Receptor,Cell-Cell Contact,all
    --species:[必需]物种信息[ human or mouse]
    --group:[可选]metadata中的组信息的slot名,默认Group
    --ident:[可选]metadata中细胞类型的slot名,默认celltype
    --outdir:[必需]输出目录
    --multiprocess:[可选]处理的线程数,默认10, [default 10] 
使用示例：
    Rscript this.r --rds rdsfile --outdir $outdir --db 'Secreted Signaling' --group Group --ident celltype --species mouse --multiprocess 10  
      \n")
  q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )      {print_usage(para) }
if (  is.null(opt$rds) )       {cat("Please give the rds file \n") ; print_usage(para)}
if ( is.null(opt$db))          {opt$db <- 'Secreted Signaling' }
if ( is.null(opt$species))     { pt$species <- 'Secreted Signaling' }
if ( is.null(opt$group))       {opt$group <- 'Group' }
if ( is.null(opt$ident))       {opt$ident <- 'celltype' }
if (  is.null(opt$outdir) )    {cat("Please give the outdir for result\n") ; print_usage(para) }
if ( is.null(opt$multiprocess)){opt$multiprocess <- 10 }

#################----------------------------------------------------------
#cellchat细胞通讯分析
#日期:2022/3/31 作者：yaojiaying
library(CellChat)
library(ggplot2)
library(patchwork) #最强大的拼图包
options(stringsAsFactors = FALSE)
library(cowplot) ###多图绘制
library(dplyr)
library(Seurat) #CombinePlots

#https://rdrr.io/github/sqjin/CellChat/man/ #函数说明

mkdirs <- function(outdir) {
    if(!file.exists(outdir)) {
        dir.create(outdir)}
    else{
        print(paste(outdir,"Dir already exists!",sep="     "))
        unlink(outdir, recursive=TRUE)
        dir.create(outdir)}
}

QC_cellchat<-function(rds,species,group.by='celltype',celltypes=names(table(rds$celltype))){
	mydata<-list(data=rds@assays$RNA@data,meta=rds@meta.data)
	mydata$meta$celltype<-as.vector(mydata$meta$celltype)
	data.input = mydata$data # normalized data matrix
	meta = mydata$meta # a dataframe with rownames containing cell mata data
	# subset部分数据进行下游分析
	#cell.use = rownames(meta)[meta$Condition == "ctrl" &　meta$Region=='PP']
	cell.use = rownames(meta)
	data.input = data.input[, cell.use]
	meta = meta[cell.use, ]
	table(meta$celltype)
	##对细胞类型排序
	meta$celltype<-factor(meta$celltype,levels=celltypes)
	cellchat <- createCellChat(object = data.input, meta = meta, group.by = group.by) #指定细胞群分群label
	return(cellchat)
}

Infer_cellchat<-function(cellchat,CellChatDB.use,thresh=0.05,thresh.p = 1,thresh.pc=0.1,min.cells = 10){
	#指定L-R database
    #预处理数据
	cellchat@DB <- CellChatDB.use
	cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
	#https://rdrr.io/github/sqjin/CellChat/man/identifyOverExpressedGenes.html
	#identifyOverExpressedGenes 不能设置成NULL,会没有var.features
	cellchat <- identifyOverExpressedGenes(cellchat,thresh.p = thresh.p,thresh.pc = thresh.pc,thresh.fc = 0,only.pos =FALSE)
	cellchat <- identifyOverExpressedInteractions(cellchat)
    #cellchat <- projectData(cellchat, PPI.human)将上一步结果保存到cellchat@LR$LRsig 中
    #计算通信概率，并推断信号网络
	cellchat <- computeCommunProb(cellchat)
	# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
	cellchat <- filterCommunication(cellchat, min.cells = min.cells)
    #计算信号通路水平上的细胞间通讯，存在net和netP中
	cellchat <- computeCommunProbPathway(cellchat,thresh = thresh)
    #聚合的细胞间通信网络
	cellchat <- aggregateNet(cellchat,thresh = thresh)
    #计算网络中心性评分
	cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
	return(cellchat)
}

Visual_cellchat_single<-function(cellchat,outdir){
	print('2.2 细胞通讯网络推断结果整理绘制circles图..')
	#整体绘图
	groupSize <- as.numeric(table(cellchat@idents))
	pdf(paste(outdir,'/netVisual_circle_all.pdf'),w=10,h=5)
	par(mfrow = c(1,2), xpd=TRUE)
	p1<-netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of Interactions")
	p2<-netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction Weights/Strength")
	print(p1)
    print(p2)
	dev.off()
	#单独可视化小图，print()不能放在一张图上
	mat <- cellchat@net$weight
	n<-ceiling(length(groupSize)/5)
	for (i in 1:nrow(mat)) {
	    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        pdf(paste(outdir,"/",mat2,'circle.pdf'),w=8,h=8)
	    mat2[i, ] <- mat[i, ]
	    p3<-netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	    print(p3) 
	}
	dev.off()
	print('2.3 单个pathway绘图（circle,hierarchy,chord,contribution,bubble图）...')
	mkdirs(paste0(outdir,'/pathways'))
	setwd(paste0(outdir,'/pathways/'))
	# Access all the signaling pathways showing significant communications
	pathways.show.all <- cellchat@netP$pathways
	# check the order of cell identity to set suitable vertex.receiver
	levels(cellchat@idents)
}

############################################ Main ############################################
outdir<-opt$outdir
prefix<-opt$prefix
rdsfile<-opt$rds
db<-opt$db
species<-opt$species
group<-opt$group
celltype<-opt$ident

setwd(outdir)
rds<-readRDS(rdsfile)

rds$Group<-rds@meta.data[,group]
rds$celltype<-rds@meta.data[,celltype]

Idents(rds)<-'celltype'
print("原始celltype数目:")
table(rds$celltype)
table(rds$Group)
groups<-names(table(rds$Group))
celltypes<-names(table(rds$celltype))
#########对细胞类型排序
rds$celltype<-factor(rds$celltype,levels=celltypes)

#1 选择数据库并且预处理表达矩阵
print('#1.选择数据库并且预处理表达矩阵...')
#（注意选择正确的物种human(1939对229种pathwayL-R互作)/mouse(2019对229种pathwayL-R互作)）
 # use CellChatDB.mouse if running on mouse data
if (species=='human'){
    CellChatDB <- CellChatDB.human
    PPI<-PPI.human
}else if (species=='mouse'){
    CellChatDB <- CellChatDB.mouse
    PPI<-PPI.mouse
}else{
    print("CellChat只能分析人和小鼠");quit()
}
dplyr::glimpse(CellChatDB$interaction)

DB<-c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
###human:1199/421/319; mouse:1209/432/378
if (db %in% DB){
    CellChatDB.use <- subsetDB(CellChatDB, search = db) # use Secreted Signaling
} else if (db =='all'){
    CellChatDB.use<-CellChatDB
} else{
    print("请选择正确的L-R互作种类：Secreted Signaling,ECM-Receptor,Cell-Cell Contact");quit()
    }

tmp_dir = paste0(outdir, '/tmp')
mkdirs(tmp_dir)
write.table(CellChatDB.use$interaction,file=paste0(outdir,'/tmp/CellChatDB.use_interaction.xls'),quote=F,sep='\t',row.names=F)
write.table(CellChatDB.human$interaction,file=paste0(outdir,'/tmp/CellChatDB.human_interaction.xls'),quote=F,sep='\t',row.names=F)
write.table(CellChatDB.mouse$interaction,file=paste0(outdir,'/tmp/CellChatDB.mouse_interaction.xls'),quote=F,sep='\t',row.names=F)

#================================================================================================
#2 按照组分别创建CellChat对象并进行细胞通讯推断分析
future::plan("multisession", workers = as.numeric(opt$multiprocess)) # do parallel

result_dir = paste0(outdir,'/1_CCI/')
mkdirs(result_dir)
print('2.按照组分别创建CellChat对象...')
#创建cellchat对象

for (i in (1:length(groups))){
    group_name = "Control"#groups[i]
    print(paste("2.1 正在处理分组 ", group_name ))
    group_outdir = paste0(result_dir,group_name)
    mkdirs(group_outdir)
    rds1<-subset(rds,Group==group_name)
    cellchat <- QC_cellchat(rds1,species,group.by='celltype',celltypes=names(table(rds1$celltype)))
    #细胞通讯推断分析
    cellchat_1a <- Infer_cellchat(cellchat,CellChatDB.use,thresh=0.05,thresh.p = 1,thresh.pc=0.1,min.cells = 10)
    net1 <- subsetCommunication(cellchat_1a,slot.name = "net",thresh = 0.05)
    netP1  <- subsetCommunication(cellchat_1a,slot.name = "netP",thresh = 0.05)
    saveRDS(cellchat_1a,file=paste0(tmp_dir,group_name,'_cellchat1.rds'))

    dataset<-CellChatDB.use$interaction
    data<-unique(dataset[,c('pathway_name','annotation')])
    netP1$annotation<- plyr::mapvalues(x =netP1$pathway_name,from = as.vector(data$pathway_name),to = as.vector(data$annotation))
    #组1结果可视化
    write.table(net1,file=paste0(group_outdir,'_net.xls'),quote=F,sep='\t',row.names=F)
    write.table(netP1,file=paste0(group_outdir,'_netP.xls'),quote=F,sep='\t',row.names=F)
    #高度可变的基因list
    var.features<-data.frame(cellchat_1a@var.features$features.info)
    var.features$gene<-rownames(var.features)
    write.table(var.features,file=paste0(group_outdir,'_celltypes_varfeatures.xls'),quote=F,sep='\t',row.names=F)

    #Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    levels(cellchat_1a@idents)
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat_1a, pattern = "outgoing")
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat_1a, pattern = "incoming")
    pdf(paste0(group_outdir,'_pathway_netP_heatmap.pdf'),w=10,h=10)
    ht1 + ht2
    dev.off()

    Visual_cellchat_single(cellchat_1a,group_outdir)

}

