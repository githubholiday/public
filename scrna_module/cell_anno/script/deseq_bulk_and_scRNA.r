library(getopt)
library(dplyr)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的rds文件',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名' ,
  'method' , 'm' , 2, 'character', '差异基因方法，默认为both ， 可选为bulk，或 scRNA , 或both',
  'cmp' , 'c', 1, 'character', '比较的组，例如："group1","group2"',
  'min.pct', 'p', 2, 'numeric', '最小百分比，默认为0.01',
    'logfc.threshold', 'l', 2, 'numeric', 'logfc阈值，默认为0.1',
    'test.use', 't', 2, 'character', '检验方法，默认为wilcox，可选为wilcox或roc等',
    'Ident' , "I" , 2 , 'character' , '聚类的标签，默认为harmony_clusters， 可选为任意列名'
),
  byrow=T,ncol=5
)

## 读取参数
args=getopt(command)

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name) || is.null(args$cmp)) {
  cat("Usage: Rscript Pseudobulk.r -i input.rds -o outdir -n name -c cmp \n")
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
print(paste0("outdir:",args$outdir))

if ( is.null(args$min.pct)) {
  min.pct <- 0.01
}else{
    min.pct <- args$min.pct
    }

print(paste0("min.pct:",min.pct))


if ( is.null(args$logfc.threshold)) {
    logfc.threshold <- 0.1
}else{
    logfc.threshold <- args$logfc.threshold
    }

print(paste0("logfc.threshold:",logfc.threshold))

if ( is.null(args$test.use)) {
    test.use <- "wilcox"
}else{
    test.use <- args$test.use
    }


print(paste0("test.use:",test.use))


if ( is.null(args$Ident)) {
    Ident.col <- "harmony_clusters"
}else{
    Ident.col <- args$Ident
}

print(paste0("Ident:",Ident.col))


output_data_frame <- function(data , output){
    data<-cbind("Gene"= rownames(data),data)
    data <- data %>%
    mutate(
      up_down = ifelse(avg_log2FC > 0, "Up", "Down"),
      significant = ifelse(p_val_adj < 0.05, "yes", "no")
    )
    write.table(data,file = output,sep = "\t",quote = FALSE,row.names = FALSE)
}



cmp <- unlist(strsplit(args$cmp , ","))
print(cmp)
if (length(cmp) != 2){
  cat("Error: cmp must be two groups!\n")
  q()
}
group1 <- cmp[1]
group2 <- cmp[2]

print(paste0("group1:",group1 , " group2:",group2))


if ( is.null(args$method)) {
  method <- "both"
}else if (args$method == "bulk"){
    method <- "bulk"
}else if (args$method == "scRNA"){
    method <- "scRNA"
}else if (args$method == "both"){
    method <- "both"
}else{
    cat("Error: method not exists!\n")
    q()
}

print(paste0("method:",method))

library(Seurat)

pbmc <- readRDS(args$input)

print(paste0("检查meta信息中是否有group和",Ident.col,"信息"))
if ( !Ident.col %in% colnames(pbmc@meta.data)){
  cat("Error: meta信息中没有",Ident.col,"信息!, 请检查\n")
  q()
}

if (! "group" %in% colnames(pbmc@meta.data)){
  cat("Error: meta信息中没有group信息!, 请检查\n")
  q()
}



if (method == "scRNA"){
  pbmc$celltype.group <- paste( unlist(pbmc[[Ident.col]]) ,  pbmc$group , sep = "_")
  Idents(pbmc) = "celltype.group"
  print(table(Idents(pbmc)))
  all_cluster <- sort(unique( unlist(pbmc[[Ident.col]]) ))
  for (i in all_cluster){
    cmp1 <- paste0(i , "_" , group1)
    cmp2 <- paste0(i , "_" , group2)
    print(paste0("cmp1:",cmp1," cmp2:",cmp2))
    print("进行scRNA差异基因分析")
    tryCatch({
      degene <- FindMarkers(pbmc, ident.1 = cmp1, ident.2 = cmp2 , min.pct = min.pct , logfc.threshold = logfc.threshold , test.use = test.use)
      prefix <- paste0(args$outdir,"/",args$name,"_scRNA_",cmp1,"_vs_",cmp2)
      print(head(degene))
      output_data_frame(degene , paste0(prefix,".xls"))
      print("完成")
    }, error = function(e){
      print(e)
      print("scRNA差异基因分析失败")
    })
  }
}else if (method == "bulk"){
    pseduo_pbmc<- AggregateExpression(pbmc, assays = "RNA", return.seurat = T, group.by = c("group","orig.ident",Ident.col))

    pseduo_pbmc$celltype.group <- paste( unlist(pseduo_pbmc[[Ident.col]]) ,  pseduo_pbmc$group , sep = "_")
    Idents(pseduo_pbmc) = "celltype.group"
  print(table(Idents(pbmc)))

    all_cluster <- sort(unique( unlist(pseduo_pbmc[[Ident.col]]) ))
    for (i in all_cluster){
      cmp1 <- paste0(i , "_" , group1)
      cmp2 <- paste0(i , "_" , group2)
      print(paste0("cmp1:",cmp1," cmp2:",cmp2))

      print("进行bulk差异基因分析")
      tryCatch({
              degene <- FindMarkers(pseduo_pbmc, ident.1 = cmp1, ident.2 = cmp2 , test.use = "DESeq2" ,min.cells.group = 1 , min.pct = min.pct , logfc.threshold = logfc.threshold)
      prefix <- paste0(args$outdir,"/",args$name,"_pesudoBulk_",cmp1,"_vs_",cmp2)
      output_data_frame(degene , paste0(prefix,".xls"))
      print("完成")

      }, error = function(e){
        print(e)
        print("bulk差异基因分析失败")
      })
    }
}else{
    pseduo_pbmc<- AggregateExpression(pbmc, assays = "RNA", return.seurat = T, group.by = c("group","orig.ident",Ident.col))
    pseduo_pbmc$celltype.group <- paste( unlist(pseduo_pbmc[[Ident.col]]) ,  pseduo_pbmc$group , sep = "_")
    Idents(pseduo_pbmc) = "celltype.group"

    pbmc$celltype.group <- paste( unlist(pbmc[[Ident.col]]) ,  pbmc$group , sep = "_")
    Idents(pbmc) = "celltype.group"
  print(table(Idents(pbmc)))

    all_cluster <- sort(unique( unlist(pseduo_pbmc[[Ident.col]]) ))



    for (i in all_cluster){

      cmp1 <- paste0(i , "_" , group1)
      cmp2 <- paste0(i , "_" , group2)
      print(paste0("cmp1:",cmp1," cmp2:",cmp2))
      print("进行bulk差异基因分析")

      cmp1_count <- table(pseduo_pbmc$celltype.group)[cmp1]
      cmp2_count <- table(pseduo_pbmc$celltype.group)[cmp2]
      if (cmp1_count > 1 && cmp2_count > 1){
        degene1 <- FindMarkers(pseduo_pbmc, ident.1 = cmp1, ident.2 = cmp2 , test.use = "DESeq2" ,min.cells.group = 1 , min.pct = min.pct , logfc.threshold = logfc.threshold)
        names(degene1) <- paste0(names(degene1), ".bulk")
        degene1$gene <- rownames(degene1)
        
        print("进行scRNA差异基因分析")
        degene2 <- FindMarkers(pbmc, ident.1 = cmp1, ident.2 = cmp2 , min.pct = min.pct , logfc.threshold = logfc.threshold , test.use = test.use)
        names(degene2) <- paste0(names(degene2), ".sc")
        degene2$gene <- rownames(degene2)

        print("合并结果")
        merge_dat <- merge(degene1 ,degene2, by = "gene")
        merge_dat <- merge_dat[order(merge_dat$p_val.sc), ]
        prefix <- paste0(args$outdir,"/",args$name,"_both_",cmp1,"_vs_",cmp2)
        output_data_frame(merge_dat , paste0( prefix, ".xls"  ) )
      }else{
        print("样本数不足，无法进行bulk 差异基因分析")
        print("只进行scRNA差异基因分析")
        tryCatch({
          degene <- FindMarkers(pbmc, ident.1 = cmp1, ident.2 = cmp2 , min.pct = min.pct , logfc.threshold = logfc.threshold , test.use = test.use)
          prefix <- paste0(args$outdir,"/",args$name,"_scRNA_",cmp1,"_vs_",cmp2)
          output_data_frame(degene , paste0(prefix,".xls"))
        }, error = function(e){
          print(e)
          print("scRNA差异基因分析失败")
        })
      }
    }
}




