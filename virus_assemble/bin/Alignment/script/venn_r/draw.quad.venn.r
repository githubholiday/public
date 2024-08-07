args = commandArgs(T)
if (length(args) != 6){
	print("Rscript draw.quad.venn.r <InFile> <OutFile><fill-color> <circle-col> <alpha>")
	print('Example : Rscript draw.quad.venn.r quad.txt quad.pdf red,blue blue,red,yellow 0.5')
	q()
	}

pdf(args[2],h=6,w=7)
data<-read.table(args[1],stringsAsFactors = FALSE,check.names = FALSE,quote = "",sep="\t")
data[,2]=as.numeric(data[,2])
library(RColorBrewer)
library(VennDiagram)
library(showtext)
showtext_begin()
font_add('SimSun', '../simsun.ttc')
k1="\n("
k2=")"
names=paste(data[1:4,1],k1,data[1:4,2],k2,sep="")
fill_col=unlist(strsplit(args[3],","))
cir_col=unlist(strsplit(args[4],","))
	venn.plot <- draw.quad.venn(
			area1 = data[1,2],
			area2 = data[2,2],
			area3 = data[3,2],
			area4 = data[4,2],
			n12 = data[5,2],
			n13 = data[6,2],
			n14 = data[7,2],
			n23 = data[8,2],
			n24 = data[9,2],
			n34 = data[10,2],
			n123 = data[11,2],
			n124 = data[12,2],
			n134 = data[13,2],
			n234 = data[14,2],
			n1234 = data[15,2],
			category = names,
			fill = fill_col,
			cex = 1.5,
			col = cir_col,
			cat.cex = 1.5,
			cat.fontfamily = "Simsun",
			cat.just = list(c(0.3,-0.1),c(0.7,0),c(0.5,0),c(0.5,0)),
			alpha=rep(args[5],4),
			lty=args[6],
			);
showtext_end()
dev.off()
