args <- commandArgs(TRUE)
if( length(args)!= 3){
    print("Rscript this.r <infile> <name> <outfile> ")
    print("Example: Rscript this.r infile.xls name outfile ")
    print("infile: 输入文件")
    print("name: 根据哪一列进行合并")
    print("outfile: 输出文件")
    q()
}


infile <- args[1]
name <- args[2]
outfile <- args[3]
library("data.table")
# X=fread('file' , sep = ' ', header = F)
dat <- read.csv(infile, sep="\t", header=T, stringsAsFactors = FALSE,check.names = FALSE)
#dat <- fread(infile, sep="\t", header=T, stringsAsFactors = FALSE,check.names = FALSE) # row.names=1, fread读取时没有这个参数。得到的对象是"data.table" "data.frame" 不能通过[,name]来取数据
# 由于第一列可能不是uniq的，所以，不能用第一列作为rowname
# dat <- read.csv(infile, sep="\t", header=T, row.names=1, stringsAsFactors = FALSE,check.names = FALSE)
dat <- as.data.frame(dat)
levels <- sort(unique(dat[,name]))
# result_dat <- matrix(1, ncol=ncol(dat)-1, nrow=length(levels))
result_dat <- matrix(1, ncol=ncol(dat), nrow=length(levels))
result_dat <- as.data.frame(result_dat)
rownames(result_dat) <- levels
# colnames(result_dat) <- colnames(dat)[2:dim(dat)[2]]
colnames(result_dat) <- colnames(dat)
for (i in levels){
	tmp_dat <- dat[dat[,name] %in% i,]
	for (j in 1:ncol(tmp_dat)){
		if (is.numeric(tmp_dat[,j])){
			result_dat[i,j] <- sum(tmp_dat[,j])
		}else{
			result_dat[i,j] <- paste(unique(tmp_dat[,j]),collapse = ",")
		}
	}
}
result_dat <- result_dat[,which(colnames(result_dat) != name) ]
result_dat <- result_dat[,which(colnames(result_dat) != 'GENE') ]
result_dat <- data.frame(name=rownames(result_dat),result_dat)
colnames(result_dat)[1] <- name
write.table(result_dat,file=outfile,sep="\t",quote=F,row.names=F)
