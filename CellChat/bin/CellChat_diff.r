library('getopt')
para<- matrix(c(
    "help"        ,"h",     0,    "logical",
    "rds"         ,"r",     2,    "character",
    "db"          ,"d",     2,    "character",
    "species"     ,"s",     2,    "character",
    "group"       ,"g",     1,    "character",
    "cmp"         ,"c",     1,    "character",
    "ident"       ,"i",     1,    "character",
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
    --cmp:[必需]提供需要分析的组名,多组时使用/分割,默认分析所有组
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
if ( is.null(opt$cmp))         {cat("Please give the cmp info \n") ; print_usage(para) }
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

DB_Select <- function( species, db='Secreted Signaling', outdir ){
    if (species=='human'){
        CellChatDB <- CellChatDB.human
        PPI<-PPI.human
    }else if (species=='mouse'){
        CellChatDB <- CellChatDB.mouse
        PPI<-PPI.mouse
    }else{
        print("CellChat只能分析人和小鼠，程序退出")
        quit()
    }
    dplyr::glimpse(CellChatDB$interaction)

    DBList<-c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
###human:1199/421/319; mouse:1209/432/378
    if (db %in% DBList){
        CellChatDB.use <- subsetDB(CellChatDB, search = db) # use Secreted Signaling
    } else if (db =='all'){
        CellChatDB.use<-CellChatDB
    } else{
        print("请选择正确的L-R互作种类：Secreted Signaling,ECM-Receptor,Cell-Cell Contact");quit()
    }

    #tmp_dir = paste0(outdir, '/tmp',sep="")
    write.table(CellChatDB.use$interaction,file=paste0(outdir,'/CellChatDB.use_interaction.xls'),quote=F,sep='\t',row.names=F)
    write.table(CellChatDB.human$interaction,file=paste0(outdir,'/CellChatDB.human_interaction.xls'),quote=F,sep='\t',row.names=F)
    write.table(CellChatDB.mouse$interaction,file=paste0(outdir,'/CellChatDB.mouse_interaction.xls'),quote=F,sep='\t',row.names=F)
    return(CellChatDB.use)
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

netVisual_pathway_plot<-function(cellchat,pathways.show,pathname,vertex.receiver=seq(1,5),sources.use = c(6:10),targets.use = c(1:5)){
#Visualize communication network associated with both signaling pathway and individual L-R pairs
	#netVisual的功能是可以同时绘制netVisual_aggregate和netVisual_individual
	#"circle", "hierarchy", "chord"三种类型图绘制只能绘制某一种
	groupSize <- as.numeric(table(cellchat@idents))
	netVisual(cellchat, signaling = pathways.show, signaling.name=pathname,vertex.receiver = vertex.receiver, layout = c("circle"), out.format ='pdf')
	netVisual(cellchat, signaling = pathways.show, signaling.name=pathname,vertex.receiver = vertex.receiver, layout = c("hierarchy"), out.format ='pdf')
	#vertex.weight = NULL 表示会按照细胞占比数绘? 不会
	# 绘制单个netVisual_individual会报错
	#netVisual(cellchat, signaling = pathways.show, signaling.name=pathname,vertex.receiver = vertex.receiver, layout = c("chord"), out.format ='pdf',vertex.weight = NULL,small.gap = 0.1,big.gap = 1)
	pdf(file =paste0(pathname,"_chord_aggregate.pdf"), width = 8, height =8)
	netVisual_aggregate(cellchat, signaling = pathways.show, signaling.name=pathname,vertex.receiver = vertex.receiver, layout = c("chord"), out.format ='pdf',vertex.weight = groupSize,small.gap = 0.1,big.gap = 1)
	dev.off()
	pdf(file =paste0(pathname,"_chord_individual.pdf"), width = 20, height =20)
	netVisual_individual(cellchat, signaling = pathways.show, signaling.name=pathname,vertex.receiver = vertex.receiver, layout = c("chord"), out.format ='pdf',vertex.weight = groupSize,small.gap = 0.1,big.gap = 1)
	dev.off()
	# Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
	#contribution barplot图绘制
	gg <- netAnalysis_contribution(cellchat, signaling = pathways.show)
	ggsave(filename=paste0(pathname, "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
	#Heatmap绘制,只能绘制单个
	p1<-netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
	pdf(paste0(pathname,'_netVisual_heatmap.pdf'),w=5,h=5)
	print(p1)
	dev.off()
	#Bubble plot 可以绘制多个绘制 #需要设定sources.use = c(6:10), targets.use = c(1:5)
	p2<-netVisual_bubble(cellchat, sources.use = sources.use, targets.use = targets.use, remove.isolate = FALSE,signaling = pathways.show)
	pdf(paste0(pathname,'_netVisual_bubble.pdf'),w=5,h=3)
	print(p2)
	dev.off()
   ##使用小提琴/点图绘制信号基因表达分布
	p1<-plotGeneExpression(cellchat, signaling = pathways.show,type='violin') #type = c("violin", "dot")
	pdf(paste0(pathname,'_plotGeneExpression_vlnplot.pdf'),w=5,h=5)
	print(p1)
	dev.off()
	#CDH,CDH2_CDH2 会有共同的基因,dotplot图会报错
	# p2<-plotGeneExpression(cellchat, signaling = pathways.show,type='dot') #type = c("violin", "dot")
	# pdf(paste0(pathname,'_plotGeneExpression_dotplot.pdf'),w=5,h=5)
	# print(p2)
	# dev.off()
}

Visual_cellchat_single<-function(cellchat,outdir){
	print('2.2 细胞通讯网络推断结果整理绘制circles图..')
	#整体绘图
	groupSize <- as.numeric(table(cellchat@idents))
	pdf(paste(outdir,'/netVisual_circle_all.pdf',sep=""),w=10,h=5)
	par(mfrow = c(1,2), xpd=TRUE)
	p1<-netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of Interactions")
	p2<-netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction Weights/Strength")
	print(p1)
    print(p2)
	dev.off()
	#单独可视化小图，print()不能放在一张图上
	mat <- cellchat@net$weight
	n<-ceiling(length(groupSize)/5)
    mat_rownames = rownames(mat)
	for (i in 1:nrow(mat)) {
	    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        pdf(paste(outdir,"/",mat_rownames[i],'circle.pdf',sep=""),w=8,h=8)
	    mat2[i, ] <- mat[i, ]
	    p3<-netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	    print(p3) 
        dev.off()
	}

	print('2.3 单个pathway绘图（circle,hierarchy,chord,contribution,bubble图）...')
	mkdirs(paste0(outdir,'/pathways'))
	setwd(paste0(outdir,'/pathways/'))
	# Access all the signaling pathways showing significant communications
	pathways.show.all <- cellchat@netP$pathways
	# check the order of cell identity to set suitable vertex.receiver
	levels(cellchat@idents)
    #levels(cellchat@idents)
    #for (i in 1:length(pathways.show.all)) {
        #pathname<-pathways.show.all[i]
        #pathways.show<-pathways.show.all[i]
        #print(paste0(i,":",pathways.show))
        #netVisual_pathway_plot(cellchat,pathways.show,pathname,vertex.receiver=vertex.receiver,sources.use=sources.use,targets.use=targets.use)
        #cellchat,pathways.show,pathname,vertex.receiver=seq(1,5),sources.use = c(6:10),targets.use = c(1:5)
	#}
}

############################################ Main ############################################
outdir<-opt$outdir
prefix<-opt$prefix
rdsfile<-opt$rds
db<-opt$db
species<-opt$species
group<-opt$group
celltype<-opt$ident
cmp <- opt$cmp
result_dir = outdir
print(paste("使用的数据库为 ",db, sep="" ))

print("读取rds文件")
rds<-readRDS(rdsfile)

rds$Group<-rds@meta.data[,group]
rds$celltype<-rds@meta.data[,celltype]

Idents(rds)<-'celltype'
print("原始celltype数目:")
table(rds$celltype)
table(rds$Group)

if(!is.null(cmp)){
    cmp_list = strsplit(cmp, "/")
    groups <- names(table(cmp_list))
}else{
    groups<-names(table(rds$Group))
}

celltypes<-names(table(rds$celltype))
#########对细胞类型排序
rds$celltype<-factor(rds$celltype,levels=celltypes)

#1 选择数据库并且预处理表达矩阵
CellChatDB.us <- DB_Select( species=species, db=db, result_dir )
#================================================================================================
#2 按照组分别创建CellChat对象并进行细胞通讯推断分析
future::plan("multisession", workers = as.numeric(opt$multiprocess)) # do parallel


mkdirs(result_dir)
print('2.按照组分别创建CellChat对象...')
#创建cellchat对象,将多个数据集合并到一起
cmp1 = groups[1]
cmp2 = groups[2]
print(paste("分析的比较组合为 ", cmp1, "*", cmp2 ))
cmp_pre = paste(cmp1, cmp2, sep="_")
result_dir = paste(outdir, "/", cmp_pre, sep="")
result_pre = paste(result_dir, "/", cmp_pre, sep="")
mkdirs(result_dir)
rds1 <- subset(rds,Group==cmp1)
rds2 <- subset(rds,Group==cmp2)
#各个样本进行cellchat
print("单样本创建cellchat对象")
cellchat1 <- QC_cellchat(rds1,species,group.by='celltype',celltypes=names(table(rds1$celltype)))
cellchat_1a <- Infer_cellchat(cellchat1,CellChatDB.use,thresh=0.05,thresh.p = 1,thresh.pc=0.1,min.cells = 10)

cellchat2 <- QC_cellchat(rds2,species,group.by='celltype',celltypes=names(table(rds1$celltype)))
cellchat_2a <- Infer_cellchat(cellchat2,CellChatDB.use,thresh=0.05,thresh.p = 1,thresh.pc=0.1,min.cells = 10)

object.list <- list(cellchat1 = cellchat1a, cellchat2 = cellchat2a)
names(object.list)<-c(cmp1,cmp2)
cellchat <- mergeCellChat(object.list, add.names = c(cmp1,cmp2))


pos.dataset<-group2
features.name<-'diff'
#组之间每个celltype做差异基因分析
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.25, thresh.p = 0.05)
net <- netMappingDEG(cellchat,features.name = features.name,thresh = 0.05)
dim(net);table(net$datasets)
str(cellchat@var.features)

net.up <- subsetCommunication(cellchat, net = net, datasets = group2)
net.down <- subsetCommunication(cellchat, net = net, datasets = group1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
write.table(net,file=paste0(result_pre,'_net_diff.xls'),quote=F,sep='\t',row.names=F)
write.table(t(net[1,]),file=paste0(result_pre,'example_net_diff.xls'),quote=F,sep='\t',row.names=T,col.names=F)

#高度可变的基因list
var.features<-data.frame(cellchat@var.features$diff.info)
var.features$gene<-rownames(var.features)
write.table(var.features,file=paste0(result_pre,'_celltypes_diffgenes.xls'),quote=F,sep='\t',row.names=F)


#存储对象
saveRDS(cellchat, file = paste0(result_pre,"_cellchat_merge.rds"))
#cellchat<-readRDS('ctrl_13Gy_cellchat_compare.rds')
saveRDS(object.list, file = paste0(result_pre,"_cellchat_object.list.rds"))


# 绘制interaction数目
library(scales)  ###颜色获取设置
color.use<-c('red','blue') #group2,group1
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(2,1),measure = 'count',color.use=color.use)
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(2,1), measure = "weight",color.use=color.use)
ggsave(filename=paste0(result_pre,'_compareInteractions.pdf'),plot=gg1 + gg2, width = 5, height = 4, units = 'in', dpi = 300)


# 红色表示在第二个组增的，蓝色表示在第一个组增加的
# ggsave和p_list一页都满足不了
pdf(paste0(result_pre,'_netVisual_diffInteraction.pdf'),w=10,h=5)
#Colors for indicating whether the signaling is increased ('color.edge[1]') or decreased ('color.edge[2]')
color.edge<-c('red','blue') #默认  color.edge = c("#b2182b", "#2166ac"),
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", color.edge = color.edge)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",color.edge = color.edge)
dev.off()

# 热图绘制
#color.heatmap A vector of two colors corresponding to max/min values, or a color name in brewer.pal only when the data in the heatmap do not contain negative values
gg1 <- netVisual_heatmap(cellchat, measure = "count",color.heatmap=c('blue','red')) #min/max values 函数说明写反了
gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.heatmap=c('blue','red'))
pdf(paste0(result_pre,'_netVisual_heatmap.pdf'),w=10,h=5)
gg1+gg2
dev.off()


# 单个interaction绘制
pdf(paste0(result_pre,'_netVisual_circle.pdf'),w=10,h=5)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

#================================================================================================
#4 定义保守和特异的pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use=c('blue','red')) ##group1,group2
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c('blue','red'))
pdf(paste0(result_pre,'_comparison_rankNet.pdf'),w=10,h=5)
gg1 + gg2
dev.off()

#5 识别上调和下调的LR pair
pdf(paste0(result_pre,'_comparison_netVisual_bubble.pdf'),w=10,h=6)
netVisual_bubble(cellchat, sources.use = sources.use,  targets.use =  targets.use,  comparison = c(1, 2), angle.x = 45,color.text=c('blue','red'))
dev.off()

pdf(paste0(result_pre,'_comparison_netVisual_bubble_DE.pdf'),w=10,h=6)
netVisual_bubble(cellchat, sources.use = sources.use,  targets.use =  targets.use, max.dataset = 1, comparison = c(1, 2), angle.x = 45,color.text=c('blue','red'))
netVisual_bubble(cellchat, sources.use = sources.use,  targets.use =  targets.use, max.dataset = 2, comparison = c(1, 2), angle.x = 45,color.text=c('blue','red'))
dev.off()

#可视化bubble
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = sources.use, targets.use = targets.use, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]),color.text=c('blue','red'))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = sources.use,  targets.use =  targets.use, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]),color.text=c('blue','red'))
#> Comparing communications on a merged object
pdf('Diff_LRpair_bubble.pdf',w=20,h=5)
gg1 + gg2
dev.off()

#可视化Chord diagram
pdf('Diff_LRpair_Chord.pdf',w=20,h=10)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = sources.use,  targets.use =  targets.use, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], sources.use = sources.use,  targets.use =  targets.use, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

print('Finish all analysis...')

