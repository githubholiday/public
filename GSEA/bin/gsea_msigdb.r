library('getopt')
para<- matrix(c(
  'help',	'h',	0,	"logical",
  'infile',	'i',	1,	"character",
  'type',	't',	2,	"character",
  'organism',	'g',	2,	"character",
  'outdir',	'o',	1,	"character"
),byrow=TRUE,ncol=4)
#===========================================================
# 仅针对人和小鼠的GO,KEGG,Reactome和MSigDB数据库，可以做enricher（传统富集方法，超几何分布，也就是Fisher exact test）和GSEA分析：
# 将基因编号转为ENTREZID(具有唯一性）：基因编号来自ANNOVAR的注释结果。
# 在转成ENTREZID时可能出现不唯一的现象。symbol与entrezID并不是绝对的一一对应的。
# 利用ClusterProfile进行富集分析：Y叔更新快，不用担心数据库过时，操作方便，出图好看易调节。。
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
  cat(getopt(para,usage=TRUE))
  cat("infile format:
      输入文件gene list(必须带有log2FC值,原则上rpkm或者order值也是可以运行成功的)，只接受SYMBOL/ENSEMBL/ENTREZID三种格式，
      暂时不接受其他格式。原则上GO数据库支持20个物种，KEGG支持几乎所有物种，Reactome支持7个物种,MSigDB支持11个物种。
      但因为org数据库包（包很大，比较占空间）只安装人和小鼠的，因此该脚本只可以做人和小鼠的基因功能富集分析，其他物种不行。
      可以做enricher（传统富集方法，超几何分布，也就是Fisher exact test）和GSEA分析。
      ===========================================================================
      *Log2FC不能为NA,lnf,-lnf值
      注意如果是lncRNA的转录组DEseq结果会有lnf,-lnf值。
      需要自行处理+1，重新得到Log2FC值,#data$Log2FC<-log2((data[,2]+1)/(data[,3]+1))
      ===========================================================================
      Gene	Log2FC	Significant
      A1BG	-4.33361161397935	yes
      A1BG-AS1	-8.23310826049641	yes
      A2M	1.03297056884149	no
      A2M-AS1	0.106991842587118	no
      A2ML1	1.0046255780959	no
      ===========================================================================
      type          enricher(ORA)和GSEA分析时候，输入文件是SYMBOL/ENSEMBL/ENTREZID。
      ===========================================================================
      outdir:outdir  of outputs,we will setwd(opt$outdir)
      Usage example:
      Rscript this.r -i gene_list_file -pr prefix -o outdir --pvalue 0.05 -t SYMBOL -g human/mouse -d GO,KEGG,Reactome,MSigDB
      Options:
      --help		h	NULL		get this help
      --infile		i	character	input gene file [forced]
      --outdir		o	character	The	resurt of out dir for analysis [forced]
      --type		t	character	one of SYMBOL/ENSEMBL/ENTREZID [forced]
      --organismg	g	character	human or mouse
      \n")
  q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )	{ print_usage(para) }
if ( is.null(opt$infile) )	{ cat("Please input the data file1 ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )	{ cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$organism) )	{ cat("Please give the organism(human/mouse) ...\n\n") ; print_usage(para) }
if ( is.null(opt$type) )	{ cat("Please give the type of gene(SYMBOL/ENSEMBL/ENTREZID) ...\n\n") ; print_usage(para) }
if ( is.null(opt$database))	{ opt$database <- 'GO,KEGG,reactome'} 
#===========================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ReactomePA)
library(msigdbr)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(ggupset)
library(cowplot)
#===========================================================
mkdirs <- function(outdir) {
  if(!file.exists(outdir)) {
    #mkdirs(dirname(fp))
    dir.create(outdir)}
  else{print(paste(outdir,"Dir already exists!",sep="     "))
    unlink(outdir, recursive=TRUE)
    dir.create(utdir,fp)
    }
}

#===========================================================
de_genes<-function(infile,species,type, outfile ){
    # infile 是基因列表文件，第一列为gene名称，需包含Log2FC列
    # species 人和小鼠的org.Hs.eg.db和org.Mm.eg.db
    # type  SYMBOL/ENSEMBL/ENTREZID
    # outfile 输出文件
    data<-read.table(infile,header=T,sep='\t')
    colnames(data)[1] <- type
    k=keys(species,keytype = type)
    geneinfo <- bitr(k,fromType = type,toType = c("ENTREZID"),OrgDb = species)
    geneset<-merge(geneinfo,data,by=type)
    #分别取差异基因构建geneList(GSEA)和gene(ORA)变量
    print(paste('enrich gene number:',length(gene)))
    #获取基因ID和log2FC列信息
    geneList <- geneset$Log2FC
    names(geneList) <- as.character(geneset$ENTREZID)
    #对其进行排序
    geneList <- sort(geneList, decreasing = TRUE)
    write.table(geneList, outfile,quote=F,sep='\t',row.names=F)
    return(geneList)
}

#===========================================================
#最多可支持"Homo sapiens,Mus musculus等11个物种，可通过msigdbr_show_species()查看。
MSigDB_clusterProfiler<-function(geneList, spe){
    m_t2g <- msigdbr(species = spe) %>% dplyr::select(gs_cat,gs_id,gs_name, entrez_gene)
    term_genes<-m_t2g[,c('gs_id','entrez_gene')]
    terms<-unique(m_t2g[,c('gs_cat','gs_id','gs_name')])
    term_names<-terms[,c('gs_id','gs_name')]
    y<-GSEA(geneList, TERM2GENE = term_genes, TERM2NAME = term_names,nPerm = 10000,pvalueCutoff =1,
    minGSSize = 1, maxGSSize = 50000)
    xy<-list( gse=y,Terms=terms)
    return(xy)
}
##只有enrich函数才可以设置readable=T获得gene symbol名
#这里统一用DOSE包的setReadable()
setReadable_write<-function(y,outfile){
    y<- setReadable(y, OrgDb = mapda,keyType = "ENTREZID")
    write.table(y,file=outfile,sep="\t", quote=FALSE, row.names=FALSE)
}

#输入数据处理
infile <- opt$infile
outdir <- opt$outdir
organism<-opt$organism
type <- opt$type

mapda<-switch(organism,
    "human"=org.Hs.eg.db,
    "mouse"=org.Mm.eg.db,
    stop("Please input the correct organism(human/mouse)!\n"))

outfile <- paste(outdir,"de_genes.txt",sep="/")
de_genes(infile,mapda,type, outfile )

opt$database <- "MSigDB"
databaseList <- strsplit(opt$database, ",")[[1]]
#5-1.MSigDB数据库(分别做ORA和GSEA test)
if ("MSigDB" %in% databaseList){
    print("Start:MSigDB analysis...")
    sign = "MSigDB"
    result_dir = paste(outdir,sign,sep="/")
    mkdirs(result_dir)

    #只能做人和小鼠的，如果不是人的，默认就是小鼠的，实际上能做的类型很多，目前不兼容
    if(organism=='human'){
        species<-'Homo sapiens'
    }else{
        species<-'"Mus musculus'
        }
    
    xyMSigDB<-MSigDB_clusterProfiler(geneList, species)
    y<- as.data.frame(setReadable(xyMSigDB$gse, OrgDb = mapda,keyType = "ENTREZID"))

    terms<-xyMSigDB$Terms[,c('gs_cat','gs_id')]
    colnames(terms)<-c('Category','ID')
    ##################################
    n2<-dim(y)[2]
    y<-merge(y, terms,by='ID')
    colN<-c('Category',colnames(y)[1:n2])
    y<-y[,colN]
    y<-y[order(y$p.adjust),]
    #为了保持sort排序保持一致，MSigDB在GSEA_plot里write，因存在大量p.adjust相同值的terms
    write.table(y,file=paste(result_dir,"gse_result.xls",sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
    saveRDS(xyMSigDB$gse, file=paste(result_dir,"gse_result.rds",sep="_"))
    print("Finish:MSigDB database with enricher(ORA) and GSEA test!")
}