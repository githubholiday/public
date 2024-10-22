library(getopt)

command=matrix(c( 
  'help',         'h', 0,'logical', '帮助文档',
  'input',        'i', 1, 'character', '输入的rds文件',
  'outdir',       'o', 1, 'character', '输出的目录',
  'name' ,        'n', 1, 'character', '输出的文件名',
  'resolution_step' , 's', 2, 'numeric', '聚类的分辨率,默认为seq(0.1, 2, 0.2)',
  'resolution' ,  'r', 2, 'numeric', '聚类的分辨率,默认为0.6',
  'method' ,      'm', 0 , 'logical', '聚类的方法,默认为 standard , 否则为SCtransform'   ,
  'integration' , 't', 2 , 'character', '整合方法，默认为Harmony ， 可选为RPCA， CCA，FastMNN',
  'nfeature' ,    'f', 2 , 'numeric', '默认选择基因的个数，默认为2000'
   ),

  byrow=T,ncol=5
  )

## 读取参数
args=getopt(command)

if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name)) {
  cat("Usage: Rscript cluster_umap -i input.rds -o outdir -n name\n")
  cat(paste(getopt(command, usage = T), "\n"))
    q() 
}

if (!file.exists(args$input)) {cat("Error: input file not exists!\n") q()}

if (!dir.exists(args$outdir)) {cat("Error: outdir not exists!\n") dir.create(args$outdir, recursive = T)}

if ( is.null(args$resolution_step)) {
  resolution_step <- seq(0.1, 2, 0.2)
}else{
   resolution_step <- as.numeric(strsplit(args$resolution_step, ",")[[1]])
}
print("聚类分辨率为"  )
print(resolution_step)

if ( is.null(args$nfeature)) {
  nfeature <- 2000
}else{
   nfeature <- as.numeric(args$nfeature)
}
print("选择的基因数为"  )
print(nfeature)

#分辨率，默认为0.6
if ( is.null(args$resolution)) {
  resolution <- 0.6
}else{
   resolution <- as.numeric(args$resolution)
}
print("聚类分辨率为"  )
print(resolution)

#标准化的方法
if ( is.null(args$method)) {
  method <- "standard"
}else{
   method <- "SCtransform"
}
print("聚类方法为"  )
print(method)

#整合方法
if ( is.null(args$integration)) {
  integration <- "Harmony"
}else{
   integration <- args$integration
}
print( paste0("整合方法为", integration)  )


#导入包
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(clustree)
library(batchelor)
library(harmony)

#输出ntop的高可变基因
output_variable_gene <- function(pbmc ,  prefix , ntop=10){
    VariableFeatures(pbmc) %>% as.data.frame %>% write.table(., file = paste(prefix,"_", ntop, "_variable_gene.txt", sep = ""), sep = "\t", quote = F, row.names = F)

    pdf(file = paste(prefix, "_variable_gene.pdf", sep = ""))
    top10 <- head(VariableFeatures(pbmc), ntop)
    plot1 <- VariableFeaturePlot(pbmc)
    print(plot1)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    print(plot2)
    dev.off()
}

output_pca_result <- function(pbmc ,  prefix){
    pdf(file = paste(prefix, "_pca.pdf", sep = ""))
    plot1 <- DimPlot(pbmc, reduction = "pca") + ggtitle("PCA plot of all samples")
    print(plot1)
    plot11 <- DimPlot(pbmc, reduction = "pca" , split.by= "orig.ident") + ggtitle("PCA plot of all samples")
    print(plot11)
    plot2 <- DimPlot(pbmc, reduction = "pca", group.by="group") + ggtitle("PCA plot of Group")
    print(plot2)
    plot2 <- DimPlot(pbmc, reduction = "pca", split.by="group") + ggtitle("PCA plot of Group")
    print(plot2)
    dev.off()

    pca_loadings <- Loadings(pbmc)
    write.csv(pca_loadings, file = paste(prefix, "_pca_loadings.csv", sep = ""), row.names = T)

    #changepoint

    pct <- pbmc[["pca"]]@stdev / sum(pbmc[["pca"]]@stdev) * 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]

    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    pcs <- min(co1, co2)
    print("拐点的pca为")
    print(pcs)

    pdf(file = paste(prefix, "_ElbowPlot.pdf", sep = ""))
    p <- ElbowPlot(pbmc , ndims = pcs+10 )#由图可看出，在PC9-PC10后出现标准差数据拐点。综上，下游分析中选择了前17个PC
    print(p)
    dev.off()

    return(pcs)    
    #https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
}

output_umap_result <- function(pbmc , prefix , method = "umap" , integrate=TRUE , seed = 0  ) {
    if (integrate){
        prefix <- paste0(prefix, "_integrate")
    }else {
        prefix <- paste0(prefix, "_unintegrate")
    }
    if (seed > 0){
        prefix <- paste0(prefix, "_seed_", seed)
    }

    pdf(file = paste0(prefix, "_" , method , ".pdf"))

    plot_name= paste0(method , "_plot")
    plot1 <- DimPlot(pbmc, reduction = method) + ggtitle(  paste0( plot_name ,  "of all samples"))
    print(plot1)
    plot11 <- DimPlot(pbmc, reduction = method , label=TRUE) + ggtitle( paste0( plot_name ,  "of all samples"))
    print(plot11)
    
    plot11 <- DimPlot(pbmc, reduction = method , split.by= "orig.ident") + ggtitle( paste0( plot_name ,  "of all samples"))
    print(plot11)
    plot2 <- DimPlot(pbmc, reduction = method, group.by="orig.ident") + ggtitle(paste0( plot_name ,  "of all samples"))
    print(plot2)
    plot2 <- DimPlot(pbmc, reduction = method, group.by="group") + ggtitle(paste0( plot_name ,  "of Group"))
    print(plot2)
    plot2 <- DimPlot(pbmc, reduction = method, split.by="group") + ggtitle(paste0( plot_name ,  "of Group"))
    print(plot2)
    dev.off()
}


do_clustree <- function(pbmc , resolution_step , prefix , integrate , method = "sct" ){
    pbmc <- FindClusters(pbmc, resolution = resolution_step )
    pdf(paste0(prefix, "_" , integrate , "_" , method, "_clustree.pdf") , width = 20 , height = 20)
    if (method == "sct") {
      p<- clustree(pbmc@meta.data, prefix = "SCT_snn_res.")
    }else{
      p<- clustree(pbmc@meta.data, prefix = "RNA_snn_res.")
    } 
    print(p)
    dev.off()
    print("绘制完clustree,请查看")
}

pbmc <- readRDS( args$input )
prefix <- paste(args$outdir, "/", args$name, sep = "")


if (method == "SCtransform") {
    ### 对线粒体比例进行回归，做SCtransform
    print("SCtransform")
    pbmc <- SCTransform(object = pbmc , vars.to.regress = "percent.mt")
    print("SCtransform完成")
    #pbmc <- RunPCA(object = pbmc , dims = 1:30)
    
    print("PCA降维")
    pbmc <- RunPCA(pbmc, 
                  ndims.print = 1:20, 
                  nfeatures.print = 50)
    best_pca <- output_pca_result(pbmc , prefix)
    choose_pca <- best_pca + 5 
    choose_pca <- max(choose_pca , 30)
    print("选择的pca为")
    print(choose_pca)
    print("PCA降维完成")

    print("计算邻接矩阵")
    pbmc <- FindNeighbors(object = pbmc, dims = 1:choose_pca)
    print("计算邻接矩阵完成")
    print("聚类")
    
    do_clustree(pbmc , resolution_step , prefix , "unintegrate" , method = "sct" )

    print( paste0("默认使用" , resolution , "的分辨率" ))
    pbmc <- FindClusters(pbmc, resolution = resolution)
    print("聚类完成")

    ### UMAP降维
    print("UMAP降维")
    pbmc <- RunUMAP(object = pbmc, dims = 1:choose_pca)
    output_umap_result(pbmc , prefix , integrate=FALSE)
    print("UMAP降维完成")

    ### TSNE降维
    print("TSNE降维")
    pbmc <- RunTSNE(object = pbmc, dims = 1:choose_pca)
    output_umap_result(pbmc , prefix , integrate=FALSE , method = "tsne")
    ### 整合

    print("整合")
    if (integration == "CCA") {
        pbmc <- IntegrateLayers(
            object = pbmc, method = CCAIntegration,
            orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method = "SCT",
            verbose = FALSE
          )
        pbmc <- RunUMAP(pbmc,  reduction="integrated.cca", reduction.name = "umap.integrated.cca" , dims = 1:choose_pca)
        pbmc <- FindNeighbors(pbmc, reduction = "integrated.cca", dims = 1:choose_pca)

        do_clustree(pbmc , resolution_step , prefix , "integrate.cca" , method = "sct" )

        pbmc <- FindClusters(pbmc, resolution = resolution, cluster.name = "cca_clusters")
        output_umap_result(pbmc , prefix , method = "umap.integrated.cca")
        
    }else if (integration == "RPCA") {
        pbmc <- IntegrateLayers(
            object = pbmc, method = RPCAIntegration,
            orig.reduction = "pca", new.reduction = "integrated.rpca",normalization.method = "SCT",
            verbose = FALSE
          )
        pbmc <- RunUMAP(pbmc,  reduction="integrated.rpca", reduction.name = "umap.integrated.rpca" , dims = 1:choose_pca)
        pbmc <- FindNeighbors(pbmc, reduction = "integrated.rpca", dims = 1:choose_pca)

        do_clustree(pbmc , resolution_step , prefix , "integrate.rpca" , method = "sct" )

        pbmc <- FindClusters(pbmc, resolution = resolution, cluster.name = "rpca_clusters")
        output_umap_result(pbmc , prefix , method = "umap.integrated.rpca")

    }else if (integration == "Harmony") {
        pbmc <- IntegrateLayers(
          object = pbmc, method = HarmonyIntegration,
          orig.reduction = "pca", new.reduction = "harmony",normalization.method = "SCT",
          verbose = FALSE
        )
        pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:choose_pca)
        do_clustree(pbmc , resolution_step , prefix , "harmony" , method = "sct" )
        pbmc <- FindClusters(pbmc, resolution = resolution, cluster.name = "harmony_clusters")
        #for ( seed in 1:100){
          #pbmc<- RunUMAP(pbmc,  reduction="harmony", reduction.name = "umap.harmony"  , dims = 1:choose_pca , seed.use= seed)
          #output_umap_result(pbmc , prefix , method = "umap.harmony" , seed=seed)
        #}
        pbmc<- RunUMAP(pbmc,  reduction="harmony", reduction.name = "umap.harmony"  , dims = 1:choose_pca )
        output_umap_result(pbmc , prefix , method = "umap.harmony" )

        pbmc <- RunTSNE(pbmc,  reduction="harmony", reduction.name = "tsne.harmony"  , dims = 1:choose_pca)
        output_umap_result(pbmc , prefix , method = "tsne.harmony")
    }else if (integration == "FastMNN") {
        pbmc <- IntegrateLayers(
          object = pbmc, method = FastMNNIntegration,
          orig.reduction = "pca", new.reduction = "integrated.mnn",normalization.method = "SCT",
          verbose = FALSE
        )
        pbmc <- RunUMAP(pbmc,  reduction="integrated.mnn", reduction.name = "umap.integrated.mnn" , dims = 1:choose_pca)
        pbmc <- FindNeighbors(pbmc, reduction = "integrated.mnn", dims = 1:choose_pca)

        do_clustree(pbmc , resolution_step , prefix , "integrated.mnn" , method = "sct" )

        pbmc <- FindClusters(pbmc, resolution = resolution, cluster.name = "mnn_clusters")
        output_umap_result(pbmc , prefix , method = "umap.integrated.mnn")

    }else{
        print("不支持的整合方法")
    }
    #DimPlot(object = pbmc, reduction = "umap")
}else{
    ### 归一化, 以10000为基数, 以log2为底
    print("归一化")
    pbmc <- NormalizeData(object = pbmc , normalization.method = "LogNormalize", scale.factor = 10000)
    print("归一化完成")
    print("鉴定高变异基因")
    ### 鉴定高变异基因，以区分细胞, 用于下游分析 , 使用VST ， 2000个基因
    pbmc <- FindVariableFeatures(object = pbmc , selection.method = "vst", nfeatures = nfeature)
    ### 查看基因的平均表达与标准方差之间的关系
    output_variable_gene(pbmc , prefix , ntop=10)
    print("鉴定高变异基因完成")

    print("标准化数据")
    ### 标准化数据 , 使用所有基因
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes)
    print("标准化数据完成")

    print("PCA降维")
    ### PCA降维 ， 使用变异基因
    pbmc <- RunPCA(pbmc, 
                  features = VariableFeatures(object = pbmc) , 
                  ndims.print = 1:20, 
                  nfeatures.print = 50)
    best_pca <- output_pca_result(pbmc , prefix)
    choose_pca <- best_pca + 5 
    choose_pca <- max(choose_pca , 20)
    print("选择的pca为")
    print(choose_pca)
    print("PCA降维完成")
    #根据JackStrawPlot来确定下一步降维作用的PC数量，如图７所示。P值越小说明PC越重要，应该纳入下一步分析
    ## 需要的cpu很多，大概30个线程
    #pbmc <- JackStraw(pbmc, num.replicate = 100)
    #pbmc <- ScoreJackStraw(pbmc, dims = 1:20)#这里dims = 1:20表示看前20个PC的P值，因为第一次算出来的结果真的惊人。而且这个dim最大只能是20.
    #JackStrawPlot(pbmc, dims = 1:20)#用JackStrawPlot绘图，绘制前10个PC的P值分布曲线
    #PCASigGenes(pbmc_small, pcs.use = 1:2)

    ### 细胞聚类
    print("计算邻接矩阵")
    pbmc <- FindNeighbors(object = pbmc, dims = 1:choose_pca)
    print("计算邻接矩阵完成")
    print("聚类")



    do_clustree(pbmc , resolution_step , prefix , "unintegrate" , method = "standard" )

    print( paste0("默认使用" , resolution , "的分辨率" ))
    pbmc <- FindClusters(pbmc, resolution = resolution)
    print("聚类完成")
    ### UMAP降维
    print("UMAP降维")
    pbmc <- RunUMAP(object = pbmc, dims = 1:choose_pca)
    output_umap_result(pbmc , prefix , integrate=FALSE)
    print("UMAP降维完成")

    ### TSNE降维
    print("tSNE降维")
    pbmc <- RunTSNE(object = pbmc, dims = 1:choose_pca)
    output_umap_result(pbmc , prefix , integrate=FALSE , method = "tsne")
    print("tSNE降维完成")

    ### 整合
    print("整合")
    if (integration == "CCA") {
        pbmc <- IntegrateLayers(
            object = pbmc, method = CCAIntegration,
            orig.reduction = "pca", new.reduction = "integrated.cca",
            verbose = FALSE
          )
        pbmc<- RunUMAP(pbmc,  reduction="integrated.cca", reduction.name = "umap.integrated.cca" , dims = 1:choose_pca)
        pbmc <- FindNeighbors(pbmc, reduction = "integrated.cca", dims = 1:choose_pca)

        do_clustree(pbmc , resolution_step , prefix , "integrate.cca" , method = "standard" )

        pbmc <- FindClusters(pbmc, resolution = resolution, cluster.name = "cca_clusters")
        output_umap_result(pbmc , prefix , method = "umap.integrated.cca")
        
    }else if (integration == "RPCA") {
        pbmc <- IntegrateLayers(
            object = pbmc, method = RPCAIntegration,
            orig.reduction = "pca", new.reduction = "integrated.rpca",
            verbose = FALSE
          )
        pbmc<- RunUMAP(pbmc,  reduction="integrated.rpca", reduction.name = "umap.integrated.rpca" , dims = 1:choose_pca)
        pbmc <- FindNeighbors(pbmc, reduction = "integrated.rpca", dims = 1:choose_pca)

        do_clustree(pbmc , resolution_step , prefix , "integrate.rpca" , method = "standard" )

        pbmc <- FindClusters(pbmc, resolution = resolution, cluster.name = "rpca_clusters")
        output_umap_result(pbmc , prefix , method = "umap.integrated.rpca")

    }else if (integration == "Harmony") {
        pbmc <- IntegrateLayers(
          object = pbmc, method = HarmonyIntegration,
          orig.reduction = "pca", new.reduction = "harmony",
          verbose = FALSE
        )
        pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:choose_pca)

        do_clustree(pbmc , resolution_step , prefix , "harmony" , method = "standard" )
        pbmc <- FindClusters(pbmc, resolution = resolution, cluster.name = "harmony_clusters")
        
        for ( seed in 1:100){
          pbmc<- RunUMAP(pbmc,  reduction="harmony", reduction.name = "umap.harmony"  , dims = 1:choose_pca , seed.use= seed)
          output_umap_result(pbmc , prefix , method = "umap.harmony" , seed=seed)
        }

        pbmc <- RunTSNE(pbmc,  reduction="harmony", reduction.name = "tsne.harmony"  , dims = 1:choose_pca)
        output_umap_result(pbmc , prefix , method = "tsne.harmony")

    }else if (integration == "FastMNN") {
        pbmc <- IntegrateLayers(
          object = pbmc, method = FastMNNIntegration,
          orig.reduction = "pca", new.reduction = "integrated.mnn",
          verbose = FALSE
        )
        pbmc<- RunUMAP(pbmc,  reduction="integrated.mnn", reduction.name = "umap.integrated.mnn" , dims = 1:choose_pca)
        pbmc <- FindNeighbors(pbmc, reduction = "integrated.mnn", dims = 1:choose_pca)

        do_clustree(pbmc , resolution_step , prefix , "integrated.mnn" , method = "standard" )

        pbmc <- FindClusters(pbmc, resolution = resolution, cluster.name = "mnn_clusters")
        output_umap_result(pbmc , prefix , method = "umap.integrated.mnn")

    }else{
        print("不支持的整合方法")
    }
    pbmc <- JoinLayers(pbmc)

}




print("保存结果")
saveRDS(pbmc, file = paste0(args$outdir, '/', args$name, '.rds'))
print("保存结果完成")



#https://mp.weixin.qq.com/s/Er-HiaeLa5usy5tdmghMlQ
# https://satijalab.org/seurat/articles/essential_commands
##https://zhuanlan.zhihu.com/p/585569143  clusttree
## https://lazappi.github.io/clustree/articles/clustree.html#using-genes-as-aesthetics

# /work/share/acuhtwkcu9/liutao/seqwisdom/9_module/scRNA/base_qc/example
# pbmc<-readRDS("outdir/rat_merge_ttclustree.rds")
# resolution =  0.6
# library(dplyr)
# library(Seurat)
# library(patchwork)
# library(ggplot2)
# pbmc <- FindClusters(pbmc, resolution = resolution)
# choose_pca <- 20
# pbmc <- RunUMAP(object = pbmc, dims = 1:choose_pca)

