library(getopt)

command=matrix(c( 
  'help', 'h', 0,'logical', '帮助文档',
  'input', 'i', 1, 'character', '输入的h5文件的配置文件，第一列为样品名，第二列为h5文件的路径',
  'outdir', 'o', 1, 'character', '输出的目录',
  'name' , 'n', 1, 'character', '输出的文件名'),
  byrow=T,ncol=5
  )
  

## 读取参数
args=getopt(command)


if (!is.null(args$help) || is.null(args$input) || is.null(args$outdir) || is.null(args$name)) {
  cat("Usage: Rscript merge_sample.r -i input.rds -o outdir -n name\n")
  cat(paste(getopt(command, usage = T), "\n"))
    q() 

}


library(dplyr)
library(Seurat)
library(patchwork)

## 读取配置文件
sample <- read.table(args$input, header = T, sep = '\t')

if ( is.null(sample$name) || is.null(sample$path) ) {
  cat("input file must have name and path column\n")
  q()
}

# 创建一个空的列表来存储Seurat对象
seurat_list <- list()

##读取每个h5文件
for (i in 1:nrow(sample)) {
    file_path <- sample$path[i]
    file_type <- sample$type[i]
  
    if (!file.exists(file_path)) {
        cat("Error: input file not exists!\n")
        print(file_path)
        q()
    }else{
        dir_info <- file.info(file_path)
    if (file_type == "10x" ) {
        print("read 10x data from dir")
        seurat_data  <- Read10X(file_path)
    }else if(file_type == "h5") {
        print("read 10x data from h5")
        seurat_data  <- Read10X_h5(file_path, use.names = TRUE, unique.features = TRUE)
    }else if(file_type == "matrix") {
        print("read 10x data from matrix")
        seurat_data  <- read.table(file_path, sep=" ")
    }else{
        print(paste("unknown file type", file_type, sep=" "))
    }
    seurat_obj   <- CreateSeuratObject(counts = seurat_data,
                                    project = sample$name[i],
                                    min.features = 200,
                                    min.cells = 3)
    seurat_list <- append(seurat_list, seurat_obj)
  }
}

## 获取所有的样品名
sample_names <- sample$name




## 合并所有的样品，添加sample.name到cell id中
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = sample_names)

### 将样品名替换成组名
name_to_group <- function(x , a_table) {
  return( a_table[ a_table$name == x, 'group' ] )
}

seurat_combined$group <- sapply(seurat_combined$orig.ident, function(x) name_to_group(x, sample) )

## 输出合并后的文件
saveRDS(seurat_combined, file = paste0(args$outdir, '/', args$name, '.rds'))


print("merge sample done!")

