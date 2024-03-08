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
library(KEGGREST)
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
  gene<-unique(as.vector(subset(geneset)$ENTREZID))
  #gene<-unique(as.vector(subset(geneset,Significant=='yes')$ENTREZID))
  print('enrich gene number:')
  print(length(gene))
  d<-geneset
  geneList <- d$Log2FC
  names(geneList) <- as.character(d$ENTREZID)
  geneList <- sort(geneList, decreasing = TRUE)
  genes<- list(gene=gene, geneList=geneList,bg=as.vector(geneinfo$ENTREZID))
  return(genes)
}
#===========================================================
#GSEA plot
#展示富集前10的条目
GSEA_plot<-function(y0,sign,database){
    library(ggpp)
    ids = c("M12056","M1519","M16362","M16843","M189","M19517","M39603","M27250")
    gsea_res <- y0
    for (k in (1:length(ids))){

        geneSetID = ids[k]
        filename<-paste(geneSetID,'GSEA_score.pdf',sep='_')
        pdf(filename, w=12, h=8)
        
        pd <- gsea_res[geneSetID, c( "NES","pvalue", "p.adjust")]
        rownames(pd) <- NULL
        for (i in seq_len(ncol(pd))){pd[, i] <- format(pd[, i], digits = 4)}
        pd_t <- data.frame(t(pd))
        data<-cbind("Type" = rownames(pd_t),pd_t)
        colnames(data) <- c("Type","Value")
        pd <- data
        p <- gseaplot2(gsea_res, geneSetID = geneSetID)

        p[[1]] <- p[[1]]+
            annotate("table", x = 12000, y = 0.4, label = pd,
            size=3,table.theme = ttheme_gtminimal
            )+
            theme(plot.title = element_text(size = 5),
                legend.position = "top",
                legend.direction = "vertical"
            )
        print(p)
        dev.off()
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
if(organism=='human'){
    mapda<-org.Hs.eg.db
}else{
    mapda<-org.Mm.eg.db
}
genes<-de_genes(indir,mapda,type)
gene<-genes$gene;geneList<-genes$geneList
bg_genes<-genes$bg

#===========================================================
#5-1.MSigDB数据库(分别做ORA和GSEA test)
if ("MSigDB" %in% databaseList){
    print("MSigDB analysis...")
    sign = "MSigDB"
    result_dir = paste(outdir,sign,sep="/")
    mkdirs(result_dir)
    setwd(result_dir)
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
    res0<-head(y,50)
    rownames(res0)<-as.vector(res0$ID)
    #5-2.ORA和GSEA结果画图
    setwd(paste(outdir,"MSigDB",sep='/'))   
    #ORA_plot(xyMSigDB$enrich,"MSigDB")
    GSEA_plot(xyMSigDB$gse,"MSigDB","MSigDB")  
    print("Finish:MSigDB database with enricher(ORA) and GSEA test!")
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

#===========================================================
