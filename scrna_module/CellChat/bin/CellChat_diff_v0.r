##### cat CellChat_diff.r 
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
  'group', 's', 2, "character",
  'cmp',   'c',    2,      "character",
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
      --group,   c,      2,      character, The  meta.data for group [default Group]
      --cmp, c, 2, character, the cmp you want to compare,A/B 
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
library(CellChat)
library(ggplot2)
library(patchwork)
options(stringsAsFactors = FALSE)
library(cowplot) ###多图绘制
library(dplyr)
library(Seurat) #CombinePlots

#https://rdrr.io/github/sqjin/CellChat/man/ #函数说明

mkdirs <- function(outdir,fp) {
  if(!file.exists(file.path(outdir,fp))) {
    #mkdirs(dirname(fp))
    dir.create(file.path(outdir,fp))}
  else{print(paste(fp,"Dir already exists!",sep="     "))
    unlink(file.path(outdir,fp), recursive=TRUE)
    dir.create(file.path(outdir,fp))}
}



QC_cellchat<-function(rds,species,group.by='celltype',celltypes=names(table(rds$celltype))){
	mydata<-list(data=rds@assays$RNA@data,meta=rds@meta.data)
	mydata$meta$celltype<-as.vector(mydata$meta$celltype)
	data.input = mydata$data # normalized data matrix
	meta = mydata$meta # a dataframe with rownames containing cell mata data
	# subset部分数据进行下游分析
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
	#CDH,CDH2_CDH2 会有共同的基因,dotplot图会报错
	# p2<-plotGeneExpression(cellchat, signaling = pathways.show,type='dot') #type = c("violin", "dot")
	# pdf(paste0(pathname,'_plotGeneExpression_dotplot.pdf'),w=5,h=5)
	# print(p2)
	# dev.off()
}



#存储所有的通路绘图结果
netVisual_pathway_plot2<-function(cellchat,pathways.show,pathname,vertex.receiver=seq(1,5),sources.use = c(6:10),targets.use = c(1:5)){
	#Visualize communication network associated with both signaling pathway and individual L-R pairs
	#netVisual的功能是可以同时绘制netVisual_aggregate和netVisual_individual
	#"circle", "hierarchy", "chord"三种类型图绘制只能绘制某一种
	groupSize<-as.numeric(table(cellchat@idents))
	#vertex.weight=groupSize #circus图
	netVisual(cellchat, signaling = pathways.show, signaling.name=pathname,vertex.receiver = vertex.receiver, layout = c("circle"), out.format ='pdf')
	netVisual(cellchat, signaling = pathways.show, signaling.name=pathname,vertex.receiver = vertex.receiver, layout = c("hierarchy"), out.format ='pdf')
	netVisual(cellchat, signaling = pathways.show, signaling.name=pathname,vertex.receiver = vertex.receiver, layout = c("chord"), out.format ='pdf')
	# Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
	#contribution barplot图绘制
	gg <- netAnalysis_contribution(cellchat, signaling = pathways.show)
	ggsave(filename=paste0(pathname, "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
	#Heatmap绘制,只能绘制单个
	# p1<-netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
	# pdf(paste0(pathname,'_netVisual_heatmap.pdf'),w=5,h=5)
	# print(p1)
	# dev.off()
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
	mkdirs(outdir,'pathways')
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



subsetCommunication1<-function(cellchat,prob=0,pval=1){
celltypes<-dimnames(cellchat@net$prob)[[1]]
LR<-dimnames(cellchat@net$prob)[[3]]
net<-data.frame()
for (i in seq(1:length(celltypes))){
	for (j in seq(1:length(celltypes))){	
			for (k in seq(1:length(LR))){	
				a<-cellchat@net$prob[i,j,k]
				b<-cellchat@net$pval[i,j,k]
				net1<-c(celltypes[i],celltypes[j],a,b,LR[k])
				net<-rbind(net,net1)
				
				}
			}
	}
colnames(net)<-c('source','target','prob','pval','interaction_name')
#source target ligand  receptor         prob pval    interaction_name
}


#data.frame(unlist(cellchat@net$prob))


outdir<-opt$outdir
prefix<-opt$prefix
rdsfile<-opt$rds
db<-opt$db
species<-opt$species
group <- opt$group
cmp<-opt$cmp
celltype<-opt$ident
#reverse<-opt$reverse

sources.use<-as.numeric(unlist(strsplit(as.character(opt$sources),split = ",",fixed=T)))
targets.use<-as.numeric(unlist(strsplit(as.character(opt$targets),split = ",",fixed=T)))
vertex.receiver<-as.numeric(unlist(strsplit(as.character(opt$vertex),split = ",",fixed=T)))
cmp_list <- strsplit(opt$cmp,split = "/",fixed=T)
print(cmp_list)

setwd(outdir)
rds<-readRDS(rdsfile)
rds$Group<-rds@meta.data[,group]
rds$celltype<-rds@meta.data[,celltype]

Idents(rds)<-'celltype'
print("原始celltype数目:")
table(rds$celltype)
table(rds$Group)
groups<-names(table(cmp_list))
print(paste0("cmp list is : ",groups[1],',',groups[2]))

if (opt$reverse =='T'){groups<-rev(groups);print(paste0("对组之间的进行reverse: ",groups))}
group1<-groups[1];group2<-groups[2]
rds1<-subset(rds,Group==group1)
rds2<-subset(rds,Group==group2)
celltypes<-names(table(rds$celltype))
#########对细胞类型排序
rds$celltype<-factor(rds$celltype,levels=celltypes)

#1 选择数据库并且预处理表达矩阵
print('#1.选择数据库并且预处理表达矩阵...')
#（注意选择正确的物种human(1939对229种pathwayL-R互作)/mouse(2019对229种pathwayL-R互作)）
 # use CellChatDB.mouse if running on mouse data
if (species=='human'){
    CellChatDB <- CellChatDB.human;PPI<-PPI.human
    }else if (species=='mouse'){
        CellChatDB <- CellChatDB.mouse;PPI<-PPI.mouse
    }else{
        print("CellChat只能分析人和小鼠");quit()}
dplyr::glimpse(CellChatDB$interaction)

DB<-c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
###human:1199/421/319; mouse:1209/432/378
if (db %in% DB){
CellChatDB.use <- subsetDB(CellChatDB, search = db) # use Secreted Signaling
} else if (db =='all'){
CellChatDB.use<-CellChatDB
} else{print("请选择正确的L-R互作种类：Secreted Signaling,ECM-Receptor,Cell-Cell Contact");quit()}


mkdirs(outdir,'tmp')
write.table(CellChatDB.use$interaction,file=paste0(outdir,'/tmp/CellChatDB.use_interaction.xls'),quote=F,sep='\t',row.names=F)
write.table(CellChatDB.human$interaction,file=paste0(outdir,'/tmp/CellChatDB.human_interaction.xls'),quote=F,sep='\t',row.names=F)
write.table(CellChatDB.mouse$interaction,file=paste0(outdir,'/tmp/CellChatDB.mouse_interaction.xls'),quote=F,sep='\t',row.names=F)

#================================================================================================
#2 按照组分别创建CellChat对象并进行细胞通讯推断分析
#future::plan("multiprocess", workers = as.numeric(opt$multiprocess)) # do parallel
future::plan("multisession", workers = as.numeric(opt$multiprocess)) # do parallel
mkdirs(outdir,'1_CCI')
setwd(paste0(outdir,'/1_CCI/'))
print('2.按照组分别创建CellChat对象...')
#创建cellchat对象
cellchat1<-QC_cellchat(rds1,species,group.by='celltype',celltypes=names(table(rds1$celltype)))
cellchat2<-QC_cellchat(rds2,species,group.by='celltype',celltypes=names(table(rds2$celltype)))
#细胞通讯推断分析
#cellchat,CellChatDB.use,thresh=0.05,thresh.p = 1,thresh.pc=0.1,min.cells = 10
cellchat1a<-Infer_cellchat(cellchat1,CellChatDB.use,thresh=0.05,thresh.p = 1,thresh.pc=0.1,min.cells = 10)
cellchat2a<-Infer_cellchat(cellchat2,CellChatDB.use,thresh=0.05,thresh.p = 1,thresh.pc=0.1,min.cells = 10)
net1<-subsetCommunication(cellchat1a,slot.name = "net",thresh = 0.05)
net2<-subsetCommunication(cellchat2a,slot.name = "net",thresh = 0.05)
netP1 <- subsetCommunication(cellchat1a,slot.name = "netP",thresh = 0.05)
netP2 <- subsetCommunication(cellchat2a,slot.name = "netP",thresh = 0.05)

saveRDS(cellchat1a,file=paste0(outdir,'/tmp/',group1,'_cellchat1.rds'))
saveRDS(cellchat2a,file=paste0(outdir,'/tmp/',group2,'_cellchat2.rds'))

dataset<-CellChatDB.use$interaction
data<-unique(dataset[,c('pathway_name','annotation')])
head(data)
netP1$annotation<- plyr::mapvalues(x =netP1$pathway_name,from = as.vector(data$pathway_name),to = as.vector(data$annotation))
netP2$annotation<- plyr::mapvalues(x =netP2$pathway_name,from = as.vector(data$pathway_name),to = as.vector(data$annotation))


#组1结果可视化
mkdirs(paste0(outdir,'/1_CCI/'),group1)
setwd(paste0(outdir,'/1_CCI/',group1))
write.table(net1,file=paste0(group1,'_net.xls'),quote=F,sep='\t',row.names=F)
write.table(netP1,file=paste0(group1,'_netP.xls'),quote=F,sep='\t',row.names=F)
#高度可变的基因list
var.features<-data.frame(cellchat1a@var.features$features.info)
var.features$gene<-rownames(var.features)
write.table(var.features,file=paste0(group1,'_celltypes_varfeatures.xls'),quote=F,sep='\t',row.names=F)

#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
n <- length(cellchat1a@netP$pathways)
levels(cellchat1a@idents)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat1a, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat1a, pattern = "incoming")
pdf(paste0(group1,'_pathway_netP_heatmap.pdf'),w=10,h=10)
ht1 + ht2
dev.off()

Visual_cellchat_single(cellchat1a,paste0(outdir,'/1_CCI/',group1),vertex.receiver=vertex.receiver,sources.use = sources.use,targets.use = targets.use)

#组2结果可视化
mkdirs(paste0(outdir,'/1_CCI/'),group2)
setwd(paste0(outdir,'/1_CCI/',group2))
write.table(net2,file=paste0(group2,'_net.xls'),quote=F,sep='\t',row.names=F)
write.table(netP2,file=paste0(group2,'_netP.xls'),quote=F,sep='\t',row.names=F)

#高度可变的基因list
var.features<-data.frame(cellchat2a@var.features$features.info)
var.features$gene<-rownames(var.features)
write.table(var.features,file=paste0(group2,'_celltypes_varfeatures.xls'),quote=F,sep='\t',row.names=F)

#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat2a, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat2a, pattern = "incoming")
pdf(paste0(group2,'_pathway_netP_heatmap.pdf'),w=10,h=10)
ht1 + ht2
dev.off()

Visual_cellchat_single(cellchat2a,paste0(outdir,'/1_CCI/',group2),vertex.receiver=vertex.receiver,sources.use = sources.use,targets.use = targets.use)

#================================================================================================
#3 组间差异比较分析
##> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
object.list <- list(cellchat1 = cellchat1a, cellchat2 = cellchat2a)
names(object.list)<-c(group1,group2)
cellchat <- mergeCellChat(object.list, add.names = c(group1,group2))
cellchat


#================================================================================================
mkdirs(outdir,'2_diff')
setwd(paste0(outdir,'/2_diff/'))
print('3.组间差异比较分析...')
#6 差异基因分析
# Identify dysfunctional signaling by using differential expression analysis
# perform differential expression analysis
#https://rdrr.io/github/sqjin/CellChat/man/identifyOverExpressedGenes.html

#identifyOverExpressedGenes
#netMappingDEG
#subsetCommunication
#extractGeneSubsetFromPair
#cellchat <- identifyOverExpressedGenes(cellchat, group.by='celltype' , only.pos = T,features.name = 'diff1', thresh.pc = 0, thresh.fc = 0, thresh.p = 1)

pos.dataset<-group2;features.name<-'diff'
#组之间每个celltype做差异基因分析
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.25, thresh.p = 0.05)
net <- netMappingDEG(cellchat,features.name = features.name,thresh = 0.05)
dim(net);table(net$datasets)
str(cellchat@var.features)

# net.up <- subsetCommunication(cellchat, net = net, datasets = group2,ligand.logFC = 0.2, receptor.logFC = NULL)
# net.down <- subsetCommunication(cellchat, net = net, datasets = group1,ligand.logFC = -0.1, receptor.logFC = -0.1)
# upregulated,downregulated
net.up <- subsetCommunication(cellchat, net = net, datasets = group2)
net.down <- subsetCommunication(cellchat, net = net, datasets = group1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
write.table(net,file=paste0(prefix,'_net_diff.xls'),quote=F,sep='\t',row.names=F)
write.table(t(net[1,]),file='example_net_diff.xls',quote=F,sep='\t',row.names=T,col.names=F)

#/annoroad/data1/bioinfo/PROJECT/big_Commercial/Cooperation/B_TET/B_TET-066/supplement/yaomengcheng/AN202106070003/Analysis/Analysis/5.celltalk/ILC_Stromal/diffgenes
#高度可变的基因list
var.features<-data.frame(cellchat@var.features$diff.info)
var.features$gene<-rownames(var.features)
write.table(var.features,file=paste0(group2,'_vs_',group1,'_celltypes_diffgenes.xls'),quote=F,sep='\t',row.names=F)


#存储对象
saveRDS(cellchat, file = paste0(outdir,'/tmp/',prefix,"_cellchat_merge.rds"))
#cellchat<-readRDS('ctrl_13Gy_cellchat_compare.rds')
saveRDS(object.list, file = paste0(outdir,'/tmp/',prefix,"_cellchat_object.list.rds"))


# 绘制interaction数目
library(scales)  ###颜色获取设置
#color.use<-hue_pal()(2)
#color.use<-c('blue','red')
color.use<-c('red','blue') #group2,group1
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(2,1),measure = 'count',color.use=color.use)
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(2,1), measure = "weight",color.use=color.use)
ggsave(filename=paste0(prefix,'_compareInteractions.pdf'),plot=gg1 + gg2, width = 5, height = 4, units = 'in', dpi = 300)


# 红色表示在第二个组增的，蓝色表示在第一个组增加的
# ggsave和p_list一页都满足不了
pdf(paste0(prefix,'_netVisual_diffInteraction.pdf'),w=10,h=5)
#Colors for indicating whether the signaling is increased ('color.edge[1]') or decreased ('color.edge[2]')
color.edge<-c('red','blue') #默认  color.edge = c("#b2182b", "#2166ac"),
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", color.edge = color.edge)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",color.edge = color.edge)
dev.off()

# 热图绘制
#color.heatmap A vector of two colors corresponding to max/min values, or a color name in brewer.pal only when the data in the heatmap do not contain negative values
gg1 <- netVisual_heatmap(cellchat, measure = "count",color.heatmap=c('blue','red')) #min/max values 函数说明写反了
gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.heatmap=c('blue','red'))
pdf(paste0(prefix,'_netVisual_heatmap.pdf'),w=10,h=5)
gg1+gg2
dev.off()


# 单个interaction绘制
pdf(paste0(prefix,'_netVisual_circle.pdf'),w=10,h=5)
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
pdf(paste0(prefix,'_comparison_rankNet.pdf'),w=10,h=5)
gg1 + gg2
dev.off()

#5 识别上调和下调的LR pair
pdf(paste0(prefix,'_comparison_netVisual_bubble.pdf'),w=10,h=6)
netVisual_bubble(cellchat, sources.use = sources.use,  targets.use =  targets.use,  comparison = c(1, 2), angle.x = 45,color.text=c('blue','red'))
dev.off()

pdf(paste0(prefix,'_comparison_netVisual_bubble_DE.pdf'),w=10,h=6)
netVisual_bubble(cellchat, sources.use = sources.use,  targets.use =  targets.use, max.dataset = 1, comparison = c(1, 2), angle.x = 45,color.text=c('blue','red'))
netVisual_bubble(cellchat, sources.use = sources.use,  targets.use =  targets.use, max.dataset = 2, comparison = c(1, 2), angle.x = 45,color.text=c('blue','red'))
dev.off()
#https://rdrr.io/github/sqjin/CellChat/man/netVisual_bubble.html
#为什么绘制不出来呢？ NA值
# LRpair<-subset(cellchat@DB$interaction,pathway_name %in% c('TGFb','IL16','IFN−II','IL2'))
# pdf(paste0(prefix,'_comparison_netVisual_bubble1.pdf'),w=10,h=6)
# netVisual_bubble(cellchat, sources.use = c(6:10), targets.use = c(1:5),  comparison = c(1, 2), angle.x = 45,color.text=c('blue','red'),pairLR.use=LRpair,thresh = 1)
# dev.off()

# Chord diagram
# pathways.show <- c("CXCL","CCL") 
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(object.list)) {
  # netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
# }

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


####单个pathway可视化绘制可视化基因小提琴图
# mkdirs(paste0(outdir,'/2_diff/'),'pathways')
# setwd(paste0(outdir,'/2_diff/pathways/'))
# cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c(group1,group2)) # set factor level
###Access all the signaling pathways showing significant communications
# pathways.show.all <- unique(c(cellchat@netP[[1]]$pathways,cellchat@netP[[2]]$pathways))
####check the order of cell identity to set suitable vertex.receiver
####vertex.receiver = seq(1,5)
####sources.use = c(6:10);targets.use = c(1:5)
# print('3 单个pathway绘图,diff（circle,hierarchy,chord,contribution,bubble图）...')
# for (i in 1:length(pathways.show.all)) {
	# pathname<-pathways.show.all[i]
	# pathways.show<-pathways.show.all[i]
	# print(paste0(i,":",pathways.show))	#netVisual_pathway_plot(cellchat,pathways.show,pathname,vertex.receiver=vertex.receiver,sources.use=sources.use,targets.use=targets.use)
	####cellchat,pathways.show,pathname,vertex.receiver=seq(1,5),sources.use = c(6:10),targets.use = c(1:5)
	# pdf(paste0(pathname,'_plotGeneExpression.pdf'),w=10,h=8)
	#如何show legend? +theme(legend.position='right')
	# p<-plotGeneExpression(cellchat, signaling = pathways.show.all[i], split.by = "datasets", colors.ggplot = T)
	# print(p)
	# dev.off()
# }

print('Finish all analysis...')
