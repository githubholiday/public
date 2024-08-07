args = commandArgs(T)
if (length(args) != 6){
	print("Rscript draw.pairwise.venn.r <InFile> <OutFile> ><fill-color> <circle-col> <alpha>")
	print("Example : Rscript draw.pairwise.venn.r pairwise.txt pairwise.pdf red,blue blue,red 0.5")
	q()
}

pdf(args[2])
d<-read.table(args[1],stringsAsFactors = FALSE,check.names = FALSE,quote = "",sep="\t")
par(font=2,font.axis=2,font.lab=2,mar=c(4,8,4,4))
library(RColorBrewer)
library(VennDiagram)
library(showtext)
showtext_begin()
font_add('SimSun', '../simsun.ttc')
k1="\n("
k2=")"
names=paste(d[1:2,1],k1,d[1:2,2],k2,sep="")
fill_col=unlist(strsplit(args[3],","))
cir_col=unlist(strsplit(args[4],","))
venn.plot <- draw.pairwise.venn(
		area1 = d[1,2],
		area2 = d[2,2],
		cross.area = d[3,2],
		category = names,
		euler.d = F,
		scaled = F,
		cat.pos = c(-5, 5),
		cat.dist = c(0.05, 0.05),
		cex = c(2,2,2),
		fill = fill_col,
		col = cir_col,
		cat.cex = c(2,2),
		alpha = rep(args[5],2),
		lty = args[6],
		cat.fontfamily = "Simsun",
		)
showtext_end()
dev.off()
