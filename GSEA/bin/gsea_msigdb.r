library('getopt')
para<- matrix(c(
  'help',	'h',	0,	"logical",
  'infile',	'i',	1,	"character",
  'prefix',	'pr',	1,	"character",
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
    cat("
【功能说明】
    该脚本进行人和小鼠的Msigdb的GSEA分析，只输出GSEA表格和rds结果，不进行绘图。如果坐其他物种，需要下载其他物种的org数据库包。
【参数说明】
    infile: 输入文件，第一列必须为gene,必须包含Log2FC列。说明:gene列的必须是SYMBOL/ENSEMBL/ENTREZID三种格式。Log2FC不能为NA,lnf,-lnf值
    prefix:输出文件前缀
    type:输入文件中gene的类型，可选 SYMBOL/ENSEMBL/ENTREZID。
    outdir:输出目录
    prefix:输出结果前缀
【使用示例】
    Usage example:
    Rscript gsea_msigdb.r -i gene_list_file --prefixr prefix -o outdir -t SYMBOL -g human/mouse
    Options:
    --help		h	NULL		get this help
    --infile		i	character	input gene file [forced]
    --outdir		o	character	The	resurt of out dir for analysis [forced]
    --prefix		pr	character	prefix of output [forced]
    --type		t	character	one of SYMBOL/ENSEMBL/ENTREZID [forced]
    --organismg	g	character	human or mouse
    \n")
    q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )	{ print_usage(para) }
if ( is.null(opt$infile) )	{ cat("Please input the data file1 ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )	{ cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$prefix) )	{ cat("Please give the prefix for outputfiles ...\n\n") ; print_usage(para) }
if ( is.null(opt$organism) )	{ cat("Please give the organism(human/mouse) ...\n\n") ; print_usage(para) }
if ( is.null(opt$type) )	{ cat("Please give the type of gene(SYMBOL/ENSEMBL/ENTREZID) ...\n\n") ; print_usage(para) }
#===========================================================
library(clusterProfiler)
library(ReactomePA)
library(msigdbr)
library(DOSE)
library(ggupset)

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
#===========================================================
mkdirs <- function(outdir) {
    if(!file.exists(outdir)) {
        dir.create(outdir)
    }else{
        print(paste(outdir,"Dir already exists!",sep="     "))
        unlink(outdir, recursive=TRUE)
        dir.create(outdir)
    }
}
#===========================================================
de_genes<-function(indir,species,type){
    #人和小鼠转化使用的物种信息
    #mapda<-org.Hs.eg.db
    #mapda<-org.Mm.eg.db
    data<-read.table(indir,header=T,sep='\t')
    colnames(data)[1] <- type
    k=keys(species,keytype = type)
    geneinfo <- bitr(k,fromType = type,toType = c("ENTREZID"),OrgDb = species)
    geneset<-merge(geneinfo,data,by=type)
    gene<-unique(as.vector(subset(geneset)$ENTREZID))
    print(paste('enrich gene number:',length(gene),sep=""))
    geneList <- geneset$Log2FC
    names(geneList) <- as.character(d$ENTREZID)
    geneList <- sort(geneList, decreasing = TRUE)
    return(geneList)
}

########################### Main ###################################

print("分析开始")
infile<-opt$infile
outdir<-opt$outdir
prefix<-opt$prefix
type<-opt$type #数据类型
organism<-opt$organism #物种，human或者mouse
#创建输出目录
mkdirs(outdir)

if(organism=='human'){
    library(org.Hs.eg.db)
    mapda<-org.Hs.eg.db
}else{
    library(org.Mm.eg.db)
    mapda<-org.Mm.eg.db
}

#1.生成geneList(for GSEA)变量
print("Start: Msigdb database with GSEA analysis")
geneList<-de_genes(infile,mapda,type)

print(head(geneList))
#只能做人和小鼠的，如果不是人的，默认就是小鼠的，实际上能做的类型很多，目前不兼容
if(organism=='human'){
    species<-'Homo sapiens'
}else{
    species<-'Mus musculus'
}

xyMSigDB<-MSigDB_clusterProfiler( geneList, species )
#将gene id转化为gene_name
y<- as.data.frame(setReadable(xyMSigDB$gse, OrgDb = mapda,keyType = "ENTREZID"))

terms<-xyMSigDB$Terms[,c('gs_cat','gs_id')]
colnames(terms)<-c('Category','ID')
n2<-dim(y)[2]
y<-merge(y, terms,by='ID')
colN<-c('Category',colnames(y)[1:n2])
y<-y[,colN]
y<-y[order(y$p.adjust),]
#为了保持sort排序保持一致，MSigDB在GSEA_plot里write，因存在大量p.adjust相同值的terms
write.table(y,file=paste(outdir,"gsea_result.xls",sep="_"),sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(xyMSigDB$gse, file=paste(outdir,"gsea_result.rds",sep="_"))

print("Finish: Msigdb database with GSEA analysis")