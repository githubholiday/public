args = commandArgs(T)
if (length(args) != 6){
	print("Rscript draw.triple.venn.r <InFile> <OutFile><fill-color> <circle-col> <alpha>")
	print("Example : Rscript draw.triple.venn.r triple.txt triple.pdf red,blue blue,red,yellow 0.5")
	q()
}

pdf(args[2],h=6,w=6)
data<-read.table(args[1],stringsAsFactors = FALSE,check.names = FALSE,quote = "",sep="\t")
data[,2]=as.numeric(data[,2])
k1=" ("
k2=")"
names=paste(data[1:3,1],k1,data[1:3,2],k2,sep="")
fill_col=unlist(strsplit(args[3],","))
cir_col=unlist(strsplit(args[4],","))

library(VennDiagram)
library(showtext)
showtext_begin()
font_add('SimSun', '../simsun.ttc')

venn.plot <- draw.triple.venn(
		area1 = data[1,2],
		area2 = data[2,2],
		area3 = data[3,2],
		n12 = data[4,2],
		n13 = data[5,2],
		n23 = data[6,2],
		n123 = data[7,2],
		category = names,
		fill = fill_col,
		col = cir_col,
		euler.d = F,
		scaled = F,
		cat.pos = c(-15,15, 180),
		cat.dist = c(0.05, 0.05, 0.025),
		cat.fontfamily = "Simsun",
		cex = 1.5,
		alpha=rep(args[5],3),
		cat.cex = 1.5,
		lty=args[6],
		)
showtext_end()
dev.off()
