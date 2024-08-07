args = commandArgs(T)
	if (length(args) != 6){
		print("Rscript draw.quintuple.venn.r <InFile> <OutFile> <fill-color> <circle-col> <alpha>")
		print("Example : Rscript draw.quintuple.venn.r quintuple.txt quintuple.pdf red,blue red.blue,red 0.5")
		q()
	}

pdf(args[2],h=9,w=9)
data<-read.table(args[1],stringsAsFactors = FALSE,check.names = FALSE,quote = "",sep="\t")
library(RColorBrewer)
library(VennDiagram)
library(showtext)
showtext_begin()
font_add('SimSun', '../simsun.ttc')
data[,2]=as.numeric(data[,2])
k1="\n("
k2=")"
names=paste(data[1:5,1],k1,data[1:5,2],k2,sep="")
fill_col=unlist(strsplit(args[3],","))
cir_col=unlist(strsplit(args[4],","))

venn.plot <- draw.quintuple.venn(
		area1 = data[1,2],
		area2 = data[2,2],
		area3 = data[3,2],
		area4 = data[4,2],
		area5 = data[5,2],
		n12 = data[6,2],
		n13 = data[7,2],
		n14 = data[8,2],
		n15 = data[9,2],
		n23 = data[10,2],
		n24 = data[11,2],
		n25 = data[12,2],
		n34 = data[13,2],
		n35 = data[14,2],
		n45 = data[15,2],
		n123 = data[16,2],
		n124 = data[17,2],
		n125 = data[18,2],
		n134 = data[19,2],
		n135 = data[20,2],
		n145 = data[21,2],
		n234 = data[22,2],
		n235 = data[23,2],
		n245 = data[24,2],
		n345 = data[25,2],
		n1234 = data[26,2],
		n1235 = data[27,2],
		n1245 = data[28,2],
		n1345 = data[29,2],
		n2345 = data[30,2],
		n12345 = data[31,2],
		#category = data[1:5,1],
		category = names,
		#fill = brewer.pal(5,"Set2"),
		fill = fill_col,
		cat.dist = c(0.38,0.25,0.5,0.5,0.25),
		#cat.col = brewer.pal(5,"Set2"),
		cat.cex = 1.4,
		cat.fontfamily = "Simsun",
		col = cir_col,
		alpha = rep(args[5], 5),
		cat.just = list(c(0.5,2),c(0.5,-0.8),c(0.5,-1.7),c(0.5,-1.5),c(0.5,-0.3)),
		margin = 0.1,
		cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.79, 1, 0.79, 1, 0.79, 1, 0.79, 1, 0.79, 1, 1, 1, 1, 1, 1.8),
		ind = TRUE,
		lty=args[6],
		);
showtext_end()
dev.off()
