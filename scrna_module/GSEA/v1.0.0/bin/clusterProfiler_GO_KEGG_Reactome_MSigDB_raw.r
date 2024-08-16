library('getopt')
para<- matrix(c(
  'help',	'h',	0,	"logical",
  'indir',	'i',	1,	"character",
  'prefix',	'pr',	1,	"character",
  'pvalue',	'p',	2,	"numeric",
  'type',	't',	2,	"character",
  'organism',	'g',	2,	"character",
  'outdir',	'o',	1,	"character",
  'database',	'd',	1,	"character"
),byrow=TRUE,ncol=4)
#===========================================================
# 仅针对人和小鼠的GO,KEGG,Reactome和MSigDB数据库，可以做enricher（传统富集方法，超几何分布，也就是Fisher exact test）和GSEA分析：
# 将基因编号转为ENTREZID(具有唯一性）：基因编号来自ANNOVAR的注释结果。
# 在转成ENTREZID时可能出现不唯一的现象。symbol与entrezID并不是绝对的一一对应的。
# 利用ClusterProfile进行富集分析：Y叔更新快，不用担心数据库过时，操作方便，出图好看易调节。。
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
  cat(getopt(para,usage=TRUE))
  cat("indir format:
      输入文件gene list(必须带有log2FC值,原则上rpkm或者order值也是可以运行成功的)，只接受SYMBOL/ENSEMBL/ENTREZID三种格式，
      暂时不接受其他格式。原则上GO数据库支持20个物种，KEGG支持几乎所有物种，Reactome支持7个物种,MSigDB支持11个物种。
      但因为org数据库包（包很大，比较占空间）只安装人和小鼠的，因此该脚本只可以做人和小鼠的基因功能富集分析，其他物种不行。
      可以做enricher（传统富集方法，超几何分布，也就是Fisher exact test）和GSEA分析。
      ===========================================================================
      *Log2FC不能为NA,lnf,-lnf值
      注意如果是lncRNA的转录组DEseq结果会有lnf,-lnf值。
      需要自行处理+1，重新得到log2FC值,#data$Log2FC<-log2((data[,2]+1)/(data[,3]+1))
      ===========================================================================
      Gene	Log2FC	Significant
      A1BG	-4.33361161397935	yes
      A1BG-AS1	-8.23310826049641	yes
      A2M	1.03297056884149	no
      A2M-AS1	0.106991842587118	no
      A2ML1	1.0046255780959	no
      ===========================================================================
      prefix        输出文件前缀
      ===========================================================================
      pvalue        GSEA画图时时候，p.adjust的阈值。
      ===========================================================================
      type          enricher(ORA)和GSEA分析时候，输入文件是SYMBOL/ENSEMBL/ENTREZID。
      ===========================================================================
      database      功能分析数据库，可选，用逗号隔开，暂时只支持 GO,KEGG,Reactome,MSigDB(如何觉得运行较慢，可以数据库分开运行)
      ===========================================================================
      outdir:outdir  of outputs,we will setwd(opt$outdir)
      Usage example:
      Rscript this.r -i gene_list_file -pr prefix -o outdir --pvalue 0.05 -t SYMBOL -g human/mouse -d GO,KEGG,Reactome,MSigDB
      Options:
      --help		h	NULL		get this help
      --indir		i	character	input gene file [forced]
      --outdir		o	character	The	resurt of out dir for analysis [forced]
      --prefix		pr	character	prefix of output [forced]
      --pvalue		pvalue	numeric	cutoff of p.adjust with GSEA plot,default=0.05
      --type		t	character	one of SYMBOL/ENSEMBL/ENTREZID [forced]
      --organismg	g	character	human or mouse
      --database	d	character	databases for ORA and GSEA test,default=GO,KEGG,reactome
      \n")
  q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )	{ print_usage(para) }
if ( is.null(opt$indir) )	{ cat("Please input the data file1 ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )	{ cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$prefix) )	{ cat("Please give the prefix for outputfiles ...\n\n") ; print_usage(para) }
if ( is.null(opt$organism) )	{ cat("Please give the organism(human/mouse) ...\n\n") ; print_usage(para) }
if ( is.null(opt$type) )	{ cat("Please give the type of gene(SYMBOL/ENSEMBL/ENTREZID) ...\n\n") ; print_usage(para) }
if ( is.null(opt$pvalue))	{ opt$pvalue <- 0.05} 
#if ( is.null(opt$type))	{ opt$type <- 'SYMBOL' }
#if ( is.null(opt$organism))	{ opt$organism <- 'human' }
if ( is.null(opt$database))	{ opt$database <- 'GO,KEGG,reactome'} 
#===========================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ReactomePA)
library(msigdbr)#MSigDB
library(DOSE)
#library(dplyr)#高效数据处理函数
library(enrichplot)
#library(gridExtra)
library(ggplot2)
library(ggupset)#upsetplot(edo)
#library(KEGGREST)#KEGGREST or reactome.db packages,KEGG.db不再用
#library(reactome.db)
#library(KEGG.db)
#library(meshes)
#library(AnnotationHub)#构建物种自己的库
#/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/R
#===========================================================
mkdirs <- function(outdir,fp) {
  if(!file.exists(file.path(outdir,fp))) {
    #mkdirs(dirname(fp))
    dir.create(file.path(outdir,fp))}
  else{print(paste(fp,"Dir already exists!",sep="     "))
    unlink(file.path(outdir,fp), recursive=TRUE)
    dir.create(file.path(outdir,fp))}
}
#===========================================================
de_genes<-function(indir,species,type){
  data<-read.table(indir,header=T,sep='\t')
  ######################
  #构建新Log2FC, 避免lnf,-lnf值
  #data$Log2FC<-log2((data[,2]+1)/(data[,3]+1))
  colnames(data)[1] <- type
  k=keys(species,keytype = type)
  #SYMBOL/ENSEMBL/ENTREZID
  geneinfo <- bitr(k,fromType = type,toType = c("ENTREZID"),OrgDb = species)
  geneset<-merge(geneinfo,data,by=type)
  #分别取差异基因构建geneList(GSEA)和gene(ORA)变量
  gene<-unique(as.vector(subset(geneset,Significant=='yes')$ENTREZID))
  print('enrich gene number:')
  print(length(gene))
  d<-geneset
  geneList <- d$Log2FC
  names(geneList) <- as.character(d$ENTREZID)
  geneList <- sort(geneList, decreasing = TRUE)
  #genes<- list(gene=gene, geneList=geneList)
  genes<- list(gene=gene, geneList=geneList,bg=as.vector(geneinfo$ENTREZID))
  return(genes)
}
#===========================================================
#GO analysis,运行较慢
GO_clusterProfiler1<-function(gene,geneList,prefix,species,ont){
    print(paste(ont,"Go terms processing..."))
    #universe= names(geneList),可设置bg gene list
    x<- enrichGO(gene= gene,OrgDb= species,pvalueCutoff  =1,qvalueCutoff  = 1,
                 minGSSize = 1, maxGSSize = 50000,ont=ont)
    y<- gseGO(geneList= geneList,OrgDb= species,
              pvalueCutoff = 1,nPerm= 10000,minGSSize= 1,maxGSSize= 50000,
              ont= ont,verbose= FALSE)
    xy<-list(enrich=x, gse=y)
  return(xy)
}
#===========================================================
#GO analysis,运行更快些
#获取level3 GO terms并对其做功能分析
#20 species,see http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
GO_clusterProfiler2<-function(gene,geneList,prefix,species,ont,nlevel,bg_genes){
  print(paste(ont,"Go terms processing..."))
  #ggo<- as.data.frame(groupGO(gene= names(geneList),OrgDb= species,ont= ont,level= 3,readable = FALSE))
  ggo<- as.data.frame(groupGO(gene= bg_genes,OrgDb= species,ont= ont,level= nlevel,readable = FALSE))
  n1<-dim(ggo)[1]
  term_list<-c()
  term_genes<-data.frame()
  for (i in (1:n1)){
    bg_gene<-unlist(strsplit(as.character(ggo$geneID[i]),split = "/",fixed=T))
    n<-ggo$Count[i]
    if (n>0){term_gene<-cbind(rep(ggo$ID[i],n),as.data.frame(bg_gene))
    colnames(term_gene)<-c("ont","gene")#rownames(term_genes)<-c(rep("GO:0030154",n))
    term_genes<-rbind(term_genes,term_gene)
    term_list<-c(term_list,as.character(ggo$ID[i]))}
  }
  term_names<-ggo[term_list,c("ID","Description")]
  x<-enricher(gene, TERM2GENE = term_genes, TERM2NAME = term_names,pvalueCutoff =1,minGSSize = 1,
              maxGSSize = 50000,qvalueCutoff=1)
  y<-GSEA(geneList, TERM2GENE = term_genes, TERM2NAME = term_names,nPerm = 10000,pvalueCutoff =1,
          minGSSize = 1, maxGSSize = 50000)
  xy<-list(enrich=x, gse=y)
  return(xy)
}
#===========================================================
#KEGG analysis
#http://www.genome.jp/kegg/catalog/org_list.html，KEGG物种默认名网站
#use_internal_data:logical, use KEGG.db(因c0008依然不能打开kegg网址，use_internal_data=T) or latest online KEGG data
KEGG_clusterProfiler<-function(gene,geneList,prefix,species,internal){
  x<-enrichKEGG(gene= gene,organism= species,pvalueCutoff =1,minGSSize = 1, maxGSSize = 50000,
                qvalueCutoff=1,use_internal_data=internal)
  y<-gseKEGG(geneList= geneList,organism= species, nPerm = 10000,minGSSize= 1,maxGSSize = 50000,
             pvalueCutoff = 1,by = "fgsea",verbose= FALSE,use_internal_data=internal)
  xy<-list(enrich=x, gse=y)
  return(xy)
}
#===========================================================
#最多可以支持human, rat, mouse, celegans, yeast,zebraﬁsh, ﬂy七个物种。
reactome_clusterProfiler<-function(gene,geneList,prefix,species){
x <- enrichPathway(gene= gene,organism= species,pvalueCutoff =1,minGSSize = 1, maxGSSize = 50000,
                   qvalueCutoff=1)
y <- gsePathway(geneList, nPerm=10000,pvalueCutoff=1,organism =species, verbose=FALSE,
                minGSSize = 1, maxGSSize = 50000)
xy<-list(enrich=x, gse=y)
return(xy)
}
#===========================================================
#最多可支持"Homo sapiens,Mus musculus等11个物种，可通过msigdbr_show_species()查看。
MSigDB_clusterProfiler<-function(gene,geneList,prefix,spe){
  m_t2g <- msigdbr(species = spe) %>% dplyr::select(gs_cat,gs_id,gs_name, entrez_gene)
  term_genes<-m_t2g[,c('gs_id','entrez_gene')]
  terms<-unique(m_t2g[,c('gs_cat','gs_id','gs_name')])
  term_names<-terms[,c('gs_id','gs_name')]
  x<-enricher(gene, TERM2GENE = term_genes, TERM2NAME = term_names,pvalueCutoff =1,minGSSize = 1,
              maxGSSize = 50000,qvalueCutoff=1)
  y<-GSEA(geneList, TERM2GENE = term_genes, TERM2NAME = term_names,nPerm = 10000,pvalueCutoff =1,
          minGSSize = 1, maxGSSize = 50000)
  xy<-list(enrich=x, gse=y,Terms=terms)
  return(xy)
}
##只有enrich函数才可以设置readable=T获得gene symbol名
#这里统一用DOSE包的setReadable()
setReadable_write<-function(x,y,sign){
  x<- setReadable(x, OrgDb =mapda,keyType = "ENTREZID")
  y<- setReadable(y, OrgDb = mapda,keyType = "ENTREZID")
  write.table(x,file=paste(prefix,sign,"enrich_result.xls",sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
  write.table(y,file=paste(prefix,sign,"gse_result.xls",sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
}

#===========================================================
#ORA plot
#展示富集前10的条目
ORA_plot<-function(x,database){
  graph <- paste(prefix,database,'enrichplot.pdf',sep='_')
  pdf(graph, w=16, h=8)
  p1<-barplot(x, showCategory=10)
  p2<-dotplot(x, showCategory=10,orderBy = "p.adjust")
  #富集到的GO terms之间的基因重叠关系进行展示
  p3<-emapplot(x,showCategory=10)
  #p4<-cnetplot(x,ategorySize="pvalue", foldChange=geneList)
  #p4<-cnetplot(x, showCategory = 10)
  p5<-upsetplot(x,showCategory=10)
  p<-plot_grid(p1,p2,p3,p5,ncol=2,labels=c("A","B","C","D"))
  print(p)
  dev.off()
}
#===========================================================
#GSEA plot
#展示富集前10的条目
GSEA_plot<-function(y0,sign,database){
  y<-y0
  graph <- paste(prefix,database,'gseplot.pdf',sep='_')
  pdf(graph, w=16, h=8)
  ID<-y$ID[1]
  main<-paste("enrichmentScore = ", signif(y[ID,]$enrichmentScore, digits = 2), sep = "")
  p1<-gseaplot2(y, geneSetID =ID,pvalue_table=T,title= main)
  p2<-ridgeplot(y,showCategory=10)
  p3<-emapplot(y,showCategory=10)
  p4<-upsetplot(y,showCategory=10)
  p<-plot_grid(p1,p2,p3,p4,ncol=2,labels=c("A","B","C","D"))
  print(p)
  dev.off()
  #提取p.adjust<0.05结果并对其画GSEA图
  cut<-pvalue
  #因y中有列变量pvalue,(p.adjust< pvalue会造成错误结果)，因此必须把pvalue替换成cut
  y<-as.data.frame(y)
  #y<-y[order(y$p.adjust),]
  #if(database=='MSigDB'){#write.table(y,file=paste(prefix,'_',sign,'_',"gse_result.xls",sep=""),sep="\t", quote=FALSE, row.names=FALSE)}
  res <- subset(y,p.adjust< cut)
  #如果显著的超过50个，只画前50个terms
  if (nrow(res) >50){
  if (database=='MSigDB'){res<-res0}else{res <- head(as.data.frame(y),50)}}
  if (nrow(res) <2){
    print ("GSEA富集的基因为0，请注意查看，将直接退出。\n")
    return
    #q()
  }
  if (nrow(res) >=2){
  print(rownames(res))
  mkdirs(outdir,paste(sign,"GSEA_plot",sep='/'))
  setwd(paste(outdir,sign,"GSEA_plot",sep='/'))
  for (k in (1:dim(res)[1])){
    i = rownames(res)[k]
    filename<-paste(i,'GSEA_score.pdf',sep='_')
    pdf(filename, w=12, h=8)
    main<-paste("enrichmentScore = ", signif(y[i,]$enrichmentScore, digits = 2), sep = "")
    p<-gseaplot2(y0, geneSetID =i,pvalue_table=T,title=main)
    # p<-gseaplot(y, geneSetID =i,title=paste(res$ID[k],"pvalue=",res$p.adjust[k]),sep=',')
    print(p)
    dev.off()}
  }
}

#===========================================================
indir<-opt$indir;prefix<-opt$prefix
pvalue<-opt$pvalue;qvalue<-opt$qvalue
type<-opt$type;organism<-opt$organism;
databaseList<-unlist(strsplit(as.character(opt$database),split = ",",fixed=T))
print(databaseList)
outdir<-opt$outdir
#mkdirs(opt$outdir,"result")
#outdir<-paste(opt$outdir,'result',sep='/')
#1.生成gene(for ORA),geneList(for GSEA)变量
print("genes build with gene(ORA(enricher)) and geneList(GSEA)...")
if(organism=='human'){mapda<-org.Hs.eg.db}else{mapda<-org.Mm.eg.db}
genes<-de_genes(indir,mapda,type)
gene<-genes$gene;geneList<-genes$geneList
bg_genes<-genes$bg
#===========================================================
#2-1.KEGG数据库(分别做ORA和GSEA test),注意internal=T(内部数据)，F(online)
#http://www.genome.jp/kegg/catalog/org_list.html，KEGG物种默认名网站
if ("KEGG" %in% databaseList){
print("KEGG analysis...")
mkdirs(outdir,"KEGG")
setwd(paste(outdir,"KEGG",sep='/'))
if(organism=='human'){species<-'hsa'}else{species<-'mmu'}
#use_internal_data:logical, use KEGG.db or latest online KEGG data
internal<-T
xyKEGG<-KEGG_clusterProfiler(gene,geneList,prefix,species,internal)
setReadable_write(xyKEGG$enrich,xyKEGG$gse,"KEGG")

#2-2.ORA和GSEA结果画图
setwd(paste(outdir,"KEGG",sep='/'))
ORA_plot(xyKEGG$enrich,"KEGG")
GSEA_plot(xyKEGG$gse,"KEGG","KEGG")
print("Finish:KEGG database with enricher(ORA) and GSEA test!")
}

#===========================================================
if ("GO" %in% databaseList){
print("GO analysis...")
#3-1.GO 数据库(分别做ORA和GSEA test)
if(organism=='human'){species<-org.Hs.eg.db}else{species<-org.Mm.eg.db}
mkdirs(outdir,"GO")
#3-2.ORA和GSEA结果画图
onts<-c("BP","CC","MF")
for(k in (1:length(onts))){
  setwd(paste(outdir,"GO",sep='/'))
  ont<-onts[k]
  mkdirs(paste(outdir,"GO",sep='/'),ont)
  setwd(paste(outdir,"GO",ont,sep='/'))
  #xyGO<-GO_clusterProfiler1(gene,geneList,prefix,species,ont) #all level
  #没有用更快的GO_clusterProfiler2是因为有时候只给50个基因，这样背景基因也是50个(已改bug)
  xyGO<-GO_clusterProfiler2(gene,geneList,prefix,species,ont,3,bg_genes) #level3
  setReadable_write(xyGO$enrich,xyGO$gse,paste("GO",ont,sep="-"))
  #xyGO<-GO_clusterProfiler3(gene,geneList,prefix,species,ont,3,bg_genes)
  #setReadable_write(xyGO$enrich,xyGO$gse,paste("GO",ont,'level3',sep="-"))
  ORA_plot(xyGO$enrich,ont)
  GSEA_plot(xyGO$gse,paste("GO",ont,sep='/'),ont)
}
print("Finish:GO database with enricher(ORA) and GSEA test!")
}

#===========================================================
#4-1.reactome数据库(分别做ORA和GSEA test)
if ("Reactome" %in% databaseList){
print("Reactome analysis...")
mkdirs(outdir,"Reactome")
setwd(paste(outdir,"Reactome",sep='/'))
if(organism=='human'){species<-'human'}else{species<-'mouse'}
xyreactome<-reactome_clusterProfiler(gene,geneList,prefix,species)
setReadable_write(xyreactome$enrich,xyreactome$gse,"Reactome")

#4-2.ORA和GSEA结果画图
setwd(paste(outdir,"Reactome",sep='/'))
ORA_plot(xyreactome$enrich,"Reactome")
GSEA_plot(xyreactome$gse,"Reactome","Reactome")
print("Finish:Reactome database with enricher(ORA) and GSEA test!")
}
#===========================================================
#5-1.MSigDB数据库(分别做ORA和GSEA test)
if ("MSigDB" %in% databaseList){
print("MSigDB analysis...")
mkdirs(outdir,"MSigDB")
setwd(paste(outdir,"MSigDB",sep='/'))
if(organism=='human'){species<-'Homo sapiens'}else{species<-'"Mus musculus'}
xyMSigDB<-MSigDB_clusterProfiler(gene,geneList,prefix,species)
#setReadable_write(xyMSigDB$enrich,xyMSigDB$gse,"MSigDB")
x<- as.data.frame(setReadable(xyMSigDB$enrich, OrgDb =mapda,keyType = "ENTREZID"))
y<- as.data.frame(setReadable(xyMSigDB$gse, OrgDb = mapda,keyType = "ENTREZID"))

terms<-xyMSigDB$Terms[,c('gs_cat','gs_id')]
colnames(terms)<-c('Category','ID')
sign<-'MSigDB'
n1<-dim(x)[2]
x<-merge(x, terms,by='ID')
colN<-c('Category',colnames(x)[1:n1])
x<-x[,colN]
##################################
n2<-dim(y)[2]
y<-merge(y, terms,by='ID')
colN<-c('Category',colnames(y)[1:n2])
y<-y[,colN]
y<-y[order(y$p.adjust),]
write.table(x,file=paste(prefix,sign,"enrich_result.xls",sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
#为了保持sort排序保持一致，MSigDB在GSEA_plot里write，因存在大量p.adjust相同值的terms
write.table(y,file=paste(prefix,sign,"gse_result.xls",sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
res0<-head(y,50)
rownames(res0)<-as.vector(res0$ID)
#5-2.ORA和GSEA结果画图
setwd(paste(outdir,"MSigDB",sep='/'))
ORA_plot(xyMSigDB$enrich,"MSigDB")
GSEA_plot(xyMSigDB$gse,"MSigDB","MSigDB")
print("Finish:MSigDB database with enricher(ORA) and GSEA test!")
}
