
#aim to plot how many genes are downregulated or upregulated upon knocking down long non-coding RNA in two different cell lines
library("VennDiagram") #loading package

######plot down-regulated genes (1 vs 0.2 and blank) 1 is NC+cisplatin; 0.2 is shRNA#1 targeting lncRNA + cisplatin; blank is shRNA#2 targeting lncRNA + cisplatin;
venn.plot <- draw.quad.venn(
  area1 = 807,
  area2 = 1109,
  area3 = 3578,
  area4 = 3669,
  n12 = 683,
  n13 = 633,
  n14 = 605,
  n23 = 776,
  n24 = 834,
  n34 = 3134,
  n123 = 561,
  n124 = 545,
  n134 = 581,
  n234 = 732,
  n1234 = 523,
  fill = c("orange", "red", "green", "blue"),
  #category = c("27(1-0.2)", "27(1-blank)", "9(1-0.2)", "9(1-blank)"), #show sample name
  
  cex=0,  # change number inside circle
  lwd=8,    #change the size of outside circle
  lty="solid"        #change the style of ourside circle, another types include dashed and blank
)



######plot up-regulated genes (1 vs 0.2 and blank)  1 is NC+cisplatin; 0.2 is shRNA#1 targeting lncRNA + cisplatin; blank is shRNA#2 targeting lncRNA + cisplatin;
venn.plot <- draw.quad.venn(
  area1 = 729,
  area2 = 938,
  area3 = 1913,
  area4 = 2108,
  n12 = 680,
  n13 = 384,
  n14 = 412,
  n23 = 454,
  n24 = 511,
  n34 = 1553,
  n123 = 363,
  n124 = 392,
  n134 = 367,
  n234 = 435,
  n1234 = 348,
  fill = c("orange", "red", "green", "blue"),
  #category = c("27(1-0.2)", "27(1-blank)", "9(1-0.2)", "9(1-blank)"), #show sample name
  
  cex=0,  # change number inside circle
  lwd=8,    #change the size of outside circle
  lty="dashed"        #change the style of ourside circle, another types include dashed and blank
)




library("pheatmap")
data=read.delim("D:/demo/song/mRNA/pathway_annotation_up.txt",head=T,sep="\t")
data2=data[,2:ncol(data)]
row.names(data2)=data[,1:1]
pheatmap(data2,border_size=1,border_color="black",show_rownames=FALSE,show_colnames=FALSE,legend=FALSE,fontsize_col=12,colorRampPalette(c("white","red"))(50),cellwidth=30,cellheight=10,cluster_rows=FALSE,cluster_cols=FALSE)


data=read.delim("D:/demo/song/mRNA/pathway_annotation_down.txt",head=T,sep="\t")
data2=data[,2:ncol(data)]
row.names(data2)=data[,1:1]
pheatmap(data2,border_size=1,border_color="black",show_rownames=FALSE,show_colnames=FALSE,legend=FALSE,fontsize_col=12,colorRampPalette(c("white","blue"))(50),cellwidth=30,cellheight=10,cluster_rows=FALSE,cluster_cols=FALSE)



