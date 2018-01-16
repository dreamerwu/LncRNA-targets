
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




library("VennDiagram")
venn.plot=draw.pairwise.venn(
  
  area1=2951,
  area2=1266,
  cross.area=38,
  #category=c("ci+","ci-"),
  cex=0,
  ext.line.lty="blank",   # design the pattern of line connecting the external area labels to their anchor points
  lwd=8,
  lty="solid",
  fill=c("red","green"),
  rotation.degree=180
)



################################
library("pheatmap")
data=read.delim("D:/demo/song/mRNA/pathway_annotation_up.txt",head=T,sep="\t")
data2=data[,2:ncol(data)]
row.names(data2)=data[,1:1]
pheatmap(data2,border_size=1,border_color="black",show_rownames=FALSE,show_colnames=FALSE,legend=FALSE,fontsize_col=12,colorRampPalette(c("white","red"))(50),cellwidth=30,cellheight=10,cluster_rows=FALSE,cluster_cols=FALSE)


data=read.delim("D:/demo/song/mRNA/pathway_annotation_down.txt",head=T,sep="\t")
data2=data[,2:ncol(data)]
row.names(data2)=data[,1:1]
pheatmap(data2,border_size=1,border_color="black",show_rownames=FALSE,show_colnames=FALSE,legend=FALSE,fontsize_col=12,colorRampPalette(c("white","blue"))(50),cellwidth=30,cellheight=10,cluster_rows=FALSE,cluster_cols=FALSE)


##############
library("ggplot2")
data=read.delim("D:/demo/song/figure/realtime.txt",head=T,sep="\t")

p=ggplot(data,aes(x=cal_27,y=scc_9,color=Gene))
pp=p+geom_point(size=5,aes(shape=group))+theme_bw()
pp+scale_color_manual(values=c("yellow","green","blue","orangered","black","maroon2","orange","gray","darkred","lightpink1",
"cyan4","darkolivegreen2","darkgreen","darkkhaki","darkgray","cyan2","cornsilk2","bisque","chocolate2","chartreuse2",
"cadetblue2","burlywood2","brown1","deepskyblue4","gold","aquamarine2","antiquewhite2","darkslategray2","darkseagreen2","darkorchid2",
"turquoise2","hotpink","thistle2","tan2","steelblue2","springgreen2","slategray2","slateblue2"))


##############January-16-2018
library("ggplot2")
ls(data)


data=read.delim("D:/demo/song/miRNA_mitochondrial/fold_change1.txt",head=T,sep="\t")
result=ggplot(data,aes(x=resistant,y=sensitive,color=group))
result2=result+geom_point(size=3)+scale_color_manual(values=c("green","gray","red"))+theme_bw()+ylim(c(-10,10))+xlim(c(-10,10))
result2+geom_abline(intercept=0,slope=1,size=1,color="blue")+theme(legend.position="none")



data=read.delim("D:/demo/song/miRNA_mitochondrial/fold_change2.txt",head=T,sep="\t")
result=ggplot(data,aes(x=total,y=mitochondrial,color=group))
result2=result+geom_point(size=3)+scale_color_manual(values=c("green","gray","red"))+theme_bw()+ylim(c(-8,8))+xlim(c(-8,8))
result2+geom_abline(intercept=0,slope=1,size=1,color="blue")+theme(legend.position="none")





data=read.delim("D:/demo/song/miRNA_mitochondrial/fold_change3.txt",head=T,sep="\t")
result=ggplot(data,aes(x=total,y=cytoplasma,color=group))
result2=result+geom_point(size=3)+scale_color_manual(values=c("green","gray","red"))+theme_bw()+ylim(c(-8,8))+xlim(c(-8,8))
result3=result2+geom_abline(intercept=0,slope=1,size=1,color="blue")+theme(legend.position="none")
result3
#result3+geom_hline(yintercept=0.585)+geom_hline(yintercept=-0.585)+geom_vline(xintercept=0.585)+geom_vline(xintercept=-0.585)


###Panle D    4 x 7.5
data=read.delim("D:/demo/song/miRNA_mitochondrial/panel_D.txt",head=T,sep="\t")
result=ggplot(data,aes(x=total_res_vs_sen,y=y))
result2=result+geom_jitter(aes(color=group_gene,size=group),height=0.1)+scale_size_manual(values=c(2,5))+scale_color_manual(values=c("orange","yellow","blue","darkblue","pink","lightblue","brown","purple","pink4","green","red","gray","black"))+theme_bw()
result2+xlim(c(-6.5,6.5))+geom_vline(xintercept=0.585,color="red",lty="dotted",lwd=2)+geom_vline(xintercept=-0.585,color="red",lty="dotted",lwd=2)+theme(legend.position="none")


###panel E 
data=read.delim("D:/demo/song/miRNA_mitochondrial/panel_E.txt",head=T,sep="\t")
result=ggplot(data,aes(x=total_res_vs_sen,y=y))
result2=result+geom_jitter(aes(color=color,size=size),height=0.1)+scale_size_manual(values=c(5,2))+scale_color_manual(values=c("steelblue2",
                                                                                                                               "yellow","green","blue","orangered","tan2","maroon2","orange","slategray2","darkred","lightpink1",
                                                                                                                               "cyan4","darkolivegreen2","darkgreen","darkkhaki","darkgray","cyan2","cornsilk2","bisque","chocolate2","chartreuse2",
                                                                                                                               "cadetblue2","burlywood2","brown1","deepskyblue4","gold","aquamarine2","antiquewhite2","darkslategray2","darkseagreen2","darkorchid2",
                                                                                                                               "turquoise2","hotpink","thistle2","gray","black"))                                                                                                                              
                                                                                                                               
                                                                                                                               
                                                                                                                               
                                                                                                                             #  "yellow","blue","darkblue","pink","lightblue","brown","purple","pink4","green","red","gray","black"))
result2+xlim(c(-6.5,6.5))+geom_vline(xintercept=0.585,color="red",lty="dotted",lwd=2)+geom_vline(xintercept=-0.585,color="red",lty="dotted",lwd=2)+theme(legend.position="none")+theme_bw()








