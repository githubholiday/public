library('getopt')
para<- matrix(c(
  'help',       'h',    0,      "logical",
  'reverse',       'r',    2,      "character",
  'rds',      'i',    1,      "character",
  'prefix',     'p',   1,      "character",
  'prob',     'b',    2,      "numeric",
  'pvalue',     'a',    2,      "numeric",
  'db',       'd',    2,      "character",
  'species',   's',    2,      "character",
  'group',   'c',    2,      "character",
  'ident',   'e',    2,      "character",
  'outdir',     'o',    1,      "character",
  'vertex',     'v',    1,      "character",
  'sources',     'l',    1,      "character",
  'targets',     't',    1,      "character",
  'multiprocess',     'm',    2,      "numeric"
),byrow=TRUE,ncol=4)
#===========================================================
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
  cat(getopt(para,usage=TRUE))
  cat("
      /annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/Rscript
      仅针对人和小鼠的物种做细胞通讯分析,注意选择正确的物种
      human(1939对229种pathwayL-R互作); mouse(2019对229种pathwayL-R互作)
      Secreted Signaling,ECM-Receptor,Cell-Cell Contact
      human:1199/421/319; mouse:1209/432/378
      1_CCI 2_diff  data
      ===========================================================================
      prefix        输出文件前缀
      ===========================================================================
      species       物种名称,human or mouse
      ===========================================================================
      db            LR互作种类（以下三种选个一个），Secreted Signaling,ECM-Receptor,Cell-Cell Contact
      ===========================================================================
      outdir:outdir  of outputs,we will setwd(opt$outdir)
      Usage example:
      Rscript this.r --rds rdsfile --outdir $outdir --prefix ILC_Stromal --db 'Secreted Signaling' --cmp Group --ident celltype --species mouse --sources 1,2,3,4,5 --targets 6,7,8,9,10 --vertex 1,2,3,4,5 --multiprocess 10 --reverse  F 
      Options:
      --help            h       NULL            get this help
      --reverse            re       charater        T(factor:condition_vs_ctrl) or F(factor:(ctrl_vs_condition))  if or not reverse cmp order (group1:ctrl)
      --rds   i       character       rds file for 10xGenomics clusters [forced]
      --outdir  o       character       The     resurt of out dir for analysis [forced]
      --prefix, pre,     1,      character,
      --pvalue, p,      2,      integer,
      --prob, pro,      2,      integer,
      --db,   d,      2,      character, [default Secreted Signaling]
      --group,   c,      2,      character, The  meta.data for cmp [default Group]
      --ident,   id,      2,      character, The  meta.data for celltype [default celltype]	  
      --species,   s,      2,      character, human or mouse
      --vertex,   v,      2,      character, order number for celltype groups 
      --sources,   l,      2,      character,order number for ligand with celltypes
      --targets,   r,      2,      character,order number for receptor with celltypes
      --multiprocess, m,      2,      integer, [default 10] 
      \n")
  q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )       { print_usage(para) }
#if ( !is.null(opt$reverse) )       { print_usage(para) }
if ( is.null(opt$rds) )       { cat("Please give the rds file ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )      { cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$prefix) )      { cat("Please give the prefix for outputfiles ...\n\n") ; print_usage(para) }
if ( is.null(opt$pvalue))       { opt$pvalue <- 0.05} 
if ( is.null(opt$prob))       { opt$prob <- 0}
if ( is.null(opt$db)) { opt$db <- 'Secreted Signaling' }
if ( is.null(opt$cmp)) { opt$cmp <- 'Group' }
if ( is.null(opt$ident)) { opt$ident <- 'celltype' }
if ( is.null(opt$species))     { cat("Please give human or mouse ...\n\n") ; print_usage(para) }
if ( is.null(opt$vertex)) { opt$sources <- '1,2,3,4,5' }
if ( is.null(opt$sources)) { opt$sources <- '1,2,3,4,5' }
if ( is.null(opt$targets)) { opt$targets <- '6,7,8,9,10' }
if ( is.null(opt$multiprocess)) { opt$multiprocess <- 10 }
if ( is.null(opt$reverse)) { opt$reverse <- 'F' }

#################----------------------------------------------------------
#cellchat细胞通讯分析
#日期:2022/3/31 作者：yaojiaying
library(CellChat)
library(ggplot2)
library(patchwork)
options(stringsAsFactors = FALSE)
library(cowplot) ###多图绘制
library(dplyr)
library(Seurat) #CombinePlots

#https://rdrr.io/github/sqjin/CellChat/man/ #函数说明

mkdirs <- function(outdir) {
  if(!file.exists(outdir)) {
    #mkdirs(dirname(fp))
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
	# cellchat<-cellchat1a
	# thresh.p = NULL;thresh.pc=NULL;min.cells = 10
	# thresh.p = 1;thresh.pc=0;min.cells = 10
	cellchat@DB <- CellChatDB.use
	cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
	#https://rdrr.io/github/sqjin/CellChat/man/identifyOverExpressedGenes.html
	#identifyOverExpressedGenes 不能设置成NULL,会没有var.features
	cellchat <- identifyOverExpressedGenes(cellchat,thresh.p = thresh.p,thresh.pc = thresh.pc,thresh.fc = 0,only.pos =FALSE)
	cellchat <- identifyOverExpressedInteractions(cellchat)
	# project gene expression data onto PPI network (optional)
	# cellchat <- projectData(cellchat, PPI)
	#source, target and ligand-receptor pair,
	cellchat <- computeCommunProb(cellchat)
	# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
	cellchat <- filterCommunication(cellchat, min.cells = min.cells)
	#df.net <- subsetCommunication(cellchat)
	#thresh = 1的cutoff是指thresh>1,而不是>=1,thresh=p值
	cellchat <- computeCommunProbPathway(cellchat,thresh = thresh)
	cellchat <- aggregateNet(cellchat,thresh = thresh)
	cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
	print(cellchat)
	print(summary(cellchat))
	print(table(cellchat@idents))
	return(cellchat)
}


#存储所有的通路绘图结果 
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

}




Visual_cellchat_single<-function(cellchat,outdir,vertex.receiver=seq(1:5),sources.use = c(6:10),targets.use = c(1:5)){
	print('2.2 细胞通讯网络推断结果整理绘图...')
	#整体绘图
	groupSize <- as.numeric(table(cellchat@idents))
	pdf('netVisual_circle_all.pdf',w=10,h=5)
	par(mfrow = c(1,2), xpd=TRUE)
	p1<-netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
	p2<-netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
	print(p1);print(p2)
	dev.off()
	#单独可视化小图，print()不能放在一张图上
	mat <- cellchat@net$weight
	n<-ceiling(length(groupSize)/5)
	pdf('netVisual_circle_part.pdf',w=8,h=8)
	#par(mfrow = c(n,5), xpd=TRUE)
	par(mfrow = c(1,1), xpd=TRUE)
	for (i in 1:nrow(mat)) {
	  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
	  mat2[i, ] <- mat[i, ]
	  p3<-netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	  print(p3) 
	}
	dev.off()
	#2.3 单个pathway绘图（"circle", "hierarchy", "chord",contribution,bubble图）
	print('2.3 单个pathway绘图（circle,hierarchy,chord,contribution,bubble图）...')
	mkdirs(paste0(outdir,'/pathways'))
	setwd(paste0(outdir,'/pathways/'))
	# Access all the signaling pathways showing significant communications
	pathways.show.all <- cellchat@netP$pathways
	# check the order of cell identity to set suitable vertex.receiver
	levels(cellchat@idents)
	# vertex.receiver = seq(1,5)
	# sources.use = c(6:10);targets.use = c(1:5)
	for (i in 1:length(pathways.show.all)) {
	  pathname<-pathways.show.all[i]
	  pathways.show<-pathways.show.all[i]
	  print(paste0(i,":",pathways.show))
	  netVisual_pathway_plot(cellchat,pathways.show,pathname,vertex.receiver=vertex.receiver,sources.use=sources.use,targets.use=targets.use)
	  #cellchat,pathways.show,pathname,vertex.receiver=seq(1,5),sources.use = c(6:10),targets.use = c(1:5)
	}
}

############################################ Main ############################################
outdir<-opt$outdir
prefix<-opt$prefix
rdsfile<-opt$rds
db<-opt$db
species<-opt$species
group<-opt$group
celltype<-opt$ident

sources.use<-as.numeric(unlist(strsplit(as.character(opt$sources),split = ",",fixed=T)))
targets.use<-as.numeric(unlist(strsplit(as.character(opt$targets),split = ",",fixed=T)))
vertex.receiver<-as.numeric(unlist(strsplit(as.character(opt$vertex),split = ",",fixed=T)))

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
if (species=='human'){CellChatDB <- CellChatDB.human;PPI<-PPI.human
}else if (species=='mouse'){CellChatDB <- CellChatDB.mouse;PPI<-PPI.mouse
}else{print("CellChat只能分析人和小鼠");quit()}
dplyr::glimpse(CellChatDB$interaction)

DB<-c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
###human:1199/421/319; mouse:1209/432/378
if (db %in% DB){
CellChatDB.use <- subsetDB(CellChatDB, search = db) # use Secreted Signaling
} else if (db =='all'){
CellChatDB.use<-CellChatDB
} else{print("请选择正确的L-R互作种类：Secreted Signaling,ECM-Receptor,Cell-Cell Contact");quit()}

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

#for (i in (1:length(groups))){
    group_name = "Control"#groups[i]
    print(paste("2.1 正在处理分组 ", group_name ))
    group_outdir = paste0(result_dir,group_name)
    mkdirs(group_outdir)
    rds1<-subset(rds,Group==group_name)
    cellchat<-QC_cellchat(rds1,species,group.by='celltype',celltypes=names(table(rds1$celltype)))
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
    n <- length(cellchat_1a@netP$pathways)
    levels(cellchat_1a@idents)
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat_1a, pattern = "outgoing")
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat_1a, pattern = "incoming")
    pdf(paste0(group_outdir,'_pathway_netP_heatmap.pdf'),w=10,h=10)
    ht1 + ht2
    dev.off()

    Visual_cellchat_single(cellchat_1a,group_outdir,vertex.receiver=vertex.receiver,sources.use = sources.use,targets.use = targets.use)

#}

