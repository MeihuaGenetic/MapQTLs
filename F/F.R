
pdf("2020-1/F.pdf",width=10,height =12)

par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,2040),ylim=c(4, 28.5+20.5+20.5+20.5+20.5-4))


#############################  I   length####################################
sub_rc <-c(100,25.5,1020, 5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

man_plot(LR_length,56,marker2_plot,marker3_plot,a=9.858573,b=11.959839)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[1],sub_rc[4]+0.5+56*20/56/4*i,sub_rc[1]+8,sub_rc[4]+0.5+56*20/56/4*i,font=2)
  text(sub_rc[1]-40,sub_rc[4]+0.5+56*20/56/4*i,i*14,cex=1.2,font=1,family="Times New Roman")
}



A<-c(1:8)
for(i in 1:8){
  segments(sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4],sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4]+0.4,font=2)
  text(sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4]-1,A[i],cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_length,9.858573),LR_length,9.858573,11.959839)$del_num
plot_text_2_length<-plot_text_2(LR_length,56,a=9.858573,marker2_plot,siggene_length_2$gene,N_del)


LR1<-plot_text_2_length[which(plot_text_2_length$Pm=="1"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:22){
  text_gene_plot2_100(LR1,i,20-0.8*i,0.6)
}

N_del<-gene_name1(signum_3(LR_length,11.959839),LR_length,9.858573,11.959839)$del_num
plot_text_3_length<-plot_text_3(LR_length,56,b=11.959839,marker3_plot,siggene_length_3$gene,N_del)


LR1<-plot_text_3_length[which(plot_text_3_length$Pm=="1"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:5){
  text_gene_plot3_250(LR1,i,8-0.8*i,0.6)
}


#text(sub_rc[1]-120,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression(-log(p-value)),srt=90,cex=1.6,family="Times New Roman")
text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-3.5,"Chromosome",cex=2,family="Times New Roman")
text(sub_rc[3]-70,sub_rc[2]-2,"I",cex=2,adj=0,family="Times New Roman")


#########################  J  daimeter  #################
sub_rc <-c(1020,25.5,1020+920, 5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

man_plot(LR_daimeter,12,marker2_plot,marker3_plot,a=7.488185 ,b=8.263014)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[3],sub_rc[4]+0.5+12*20/12/4*i,sub_rc[3]-8,sub_rc[4]+0.5+12*20/12/4*i,font=2)
  text(sub_rc[3]+40,sub_rc[4]+0.5+12*20/12/4*i,i*3,cex=1.2,font=1,family="Times New Roman")
}


A<-c(1:8)
for(i in 1:8){
  segments(sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4],sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4]+0.4,font=2)
  text(sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4]-1,A[i],cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_daimeter,7.488185),LR_daimeter,7.488185 ,8.263014)$del_num
plot_text_2_daimeter<-plot_text_2(LR_daimeter,12,a=7.488185,marker2_plot,siggene_daimeter_2$gene,N_del)

LR1<-plot_text_2_daimeter[which(plot_text_2_daimeter$Pm=="2"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:6){
  text_gene_plot2_300(LR1,i,20-0.8*i,0.6)
}


LR1<-plot_text_2_daimeter[which(plot_text_2_daimeter$Pm=="6"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:5){
  text_gene_plot2_550(LR1,i,20-0.8*i,0.6)
}

text(plot_text_2_daimeter$cm[12]+sub_rc[1]+20+20*(plot_text_2_daimeter$Pm[12]-1),
     plot_text_2_daimeter$LR[12]+sub_rc[4]+0.5,plot_text_2_daimeter$gene[12],adj=0,cex=0.6,family="Times New Roman1")



N_del<-gene_name1(signum_3(LR_daimeter,8.263014),LR_daimeter,7.488185 ,8.263014)$del_num
plot_text_3_daimeter<-plot_text_3(LR_daimeter,12,b=8.263014,marker3_plot,siggene_daimeter_3$gene,N_del)

text(plot_text_3_daimeter$cm[1]+sub_rc[1]+20+20*(plot_text_3_daimeter$Pm[1]-1),
     plot_text_3_daimeter$LR[1]+sub_rc[4]+0.5,plot_text_3_daimeter$gene[1],adj=0,cex=0.6,family="Times New Roman1")

text(plot_text_3_daimeter$cm[2]+sub_rc[1]+20+20*(plot_text_3_daimeter$Pm[2]-1),
     plot_text_3_daimeter$LR[2]+sub_rc[4]+0.5,plot_text_3_daimeter$gene[2],adj=0,cex=0.6,family="Times New Roman1")

text(plot_text_3_daimeter$cm[3]+sub_rc[1]+20+20*(plot_text_3_daimeter$Pm[3]-1),
     plot_text_3_daimeter$LR[3]+sub_rc[4]+0.5,plot_text_3_daimeter$gene[3],adj=0,cex=0.6,family="Times New Roman1")


text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-3.5,"Chromosome",cex=2,family="Times New Roman")
text(sub_rc[1]+10,sub_rc[2]-2,"J",cex=2,adj=0,family="Times New Roman")


#############################  G  colour####################################
sub_rc <-c(100,25.5+20.5,1020, 5+20.5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

man_plot(LR_colour,72,marker2_plot,marker3_plot,a=10.768008,b=7.654453)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:3){
  segments(sub_rc[1],sub_rc[4]+0.5+72*20/72/4*i,sub_rc[1]+8,sub_rc[4]+0.5+72*20/72/4*i,font=2)
  text(sub_rc[1]-40,sub_rc[4]+0.5+72*20/72/4*i,i*18,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_colour,10.768008),LR_colour,10.768008,7.654453)$del_num
plot_text_2_colour<-plot_text_2(LR_colour,72,a=10.768008,marker2_plot,siggene_colour_2$gene,N_del)


LR1<-plot_text_2_colour[which(plot_text_2_colour$Pm=="4"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]


for(i in 1:13){
  text_gene_plot2_600(LR1,i,20-0.8*i,0.6)
}


N_del<-gene_name1(signum_3(LR_colour,7.654453),LR_colour,10.768008,7.654453)$del_num
plot_text_3_colour<-plot_text_3(LR_colour,72,b=7.654453,marker3_plot,siggene_colour_3$gene,N_del)


for(i in 1:2){
  text_gene_plot3(plot_text_3_colour,i,1,0.6)
}
text_gene_plot3(plot_text_3_colour,3,0,0.6)

text(sub_rc[1]+10,sub_rc[2]-2,"G",cex=2,adj=0,family="Times New Roman")


#########################  H  period  #################
sub_rc <-c(1020,25.5+20.5,1020+920, 5+20.5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

man_plot(LR_period,10,marker2_plot,marker3_plot,a=6.878382,b=6.168139)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:3){
  segments(sub_rc[3],sub_rc[4]+0.5+9*20/10/3*i,sub_rc[3]-8,sub_rc[4]+0.5+9*20/10/3*i,font=2)
  text(sub_rc[3]+40,sub_rc[4]+0.5+9*20/10/3*i,i*3,cex=1.2,font=1,family="Times New Roman")
}


N_del<-gene_name1(signum_2(LR_period,6.878382),LR_period,6.878382,6.168139)$del_num
plot_text_2_period<-plot_text_2(LR_period,10,a=6.878382,marker2_plot,siggene_period_2$gene,N_del)


LR1<-plot_text_2_period[which(plot_text_2_period$Pm=="1"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:4){
  text_gene_plot2_100(LR1,i,20-0.8*i,0.6)
}


text(sub_rc[3]-70,sub_rc[2]-2,"H",cex=2,adj=0,family="Times New Roman")


#############################  E  pistils####################################
sub_rc <-c(100,25.5+20.5+20.5,1020, 5+20.5+20.5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

man_plot(LR_pistils,40,marker2_plot,marker3_plot,a=14.983331,b=10.665697)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[1],sub_rc[4]+0.5+40*20/40/4*i,sub_rc[1]+8,sub_rc[4]+0.5+40*20/40/4*i,font=2)
  text(sub_rc[1]-40,sub_rc[4]+0.5+40*20/40/4*i,i*10,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_pistils,14.983331),LR_pistils,14.983331,10.665697)$del_num
plot_text_2_pistils<-plot_text_2(LR_pistils,40,a=14.983331,marker2_plot,siggene_pistils_2$gene,N_del)

LR1<-plot_text_2_pistils[which(plot_text_2_pistils$Pm=="1"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:17){
  text_gene_plot2_100(LR1,i,20-0.8*i,0.6)
}

N_del<-gene_name1(signum_3(LR_pistils,10.665697),LR_pistils,14.983331,10.665697)$del_num
plot_text_3_pistils<-plot_text_3(LR_pistils,40,b=10.665697,marker3_plot,siggene_pistils_3$gene,N_del)


LR1<-plot_text_3_pistils[which(plot_text_3_pistils$Pm=="1"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:17){
  text_gene_plot3_250(LR1,i,9-0.8*i,0.6)
}

text(sub_rc[1]-120,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression(-log(p-value)),srt=90,cex=1.6,family="Times New Roman")


text(sub_rc[3]-70,sub_rc[2]-2,"E",cex=2,adj=0,family="Times New Roman")


#########################  F  bud  #################
sub_rc <-c(1020,25.5+20.5+20.5,1020+920, 5+20.5+20.5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

man_plot(LR_bud,7.4,marker2_plot,marker3_plot,a=5.634593,b=5.248507)


rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[3],sub_rc[4]+0.5+6*20/7.4/3*i,sub_rc[3]-8,sub_rc[4]+0.5+6*20/7.4/3*i,font=2)
  text(sub_rc[3]+40,sub_rc[4]+0.5+6*20/7.4/3*i,i*2,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_bud,5.634593),LR_bud,5.634593,5.248507)$del_num
plot_text_2_bud<-plot_text_2(LR_bud,a=5.634593,marker2_plot,siggene_bud_2$gene,N_del)


text(plot_text_2_bud$cm[1]+sub_rc[1]+20+20*(plot_text_2_bud$Pm[1]-1),
     plot_text_2_bud$LR[1]+sub_rc[4]+0.5,plot_text_2_bud$gene[1],adj=1,cex=0.6,family="Times New Roman1")

text(plot_text_2_bud$cm[2]+sub_rc[1]+20+20*(plot_text_2_bud$Pm[2]-1),
     plot_text_2_bud$LR[2]+sub_rc[4]+0.5,plot_text_2_bud$gene[2],adj=0,cex=0.6,family="Times New Roman1")

text(plot_text_2_bud$cm[3]+sub_rc[1]+20+20*(plot_text_2_bud$Pm[3]-1),
     plot_text_2_bud$LR[3]+sub_rc[4]+0.5,plot_text_2_bud$gene[3],adj=1,cex=0.6,family="Times New Roman1")

text(plot_text_2_bud$cm[4]+sub_rc[1]+30+20*(plot_text_2_bud$Pm[4]-1),
     plot_text_2_bud$LR[4]+sub_rc[4]+1,plot_text_2_bud$gene[4],adj=0,cex=0.6,family="Times New Roman1")

segments(plot_text_2_bud$cm[4]+sub_rc[1]+20+20*(plot_text_2_bud$Pm[4]-1),
         plot_text_2_bud$LR[4]+sub_rc[4]+0.5,
         plot_text_2_bud$cm[4]+sub_rc[1]+30+20*(plot_text_2_bud$Pm[4]-1),
         plot_text_2_bud$LR[4]+sub_rc[4]+1,lwd=0.5)


text(sub_rc[3]-70,sub_rc[2]-2,"F",cex=2,adj=0,family="Times New Roman")


text(sub_rc[3]+120,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression(-log(p-value)),srt=90,cex=1.6,family="Times New Roman")



#############################  C  shape####################################
sub_rc <-c(100,25.5+20.5+20.5+20.5,1020, 5+20.5+20.5+20.5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

man_plot(LR_shape,9.2,marker2_plot,marker3_plot,a=7.040742,b=5.254234)


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+0.5+8*20/9.2/4*i,sub_rc[1]+8,sub_rc[4]+0.5+8*20/9.2/4*i,font=2)
  text(sub_rc[1]-40,sub_rc[4]+0.5+8*20/9.2/4*i,i*2,cex=1.2,font=1,family="Times New Roman")
}


N_del<-gene_name1(signum_2(LR_shape,7.040742),LR_shape,7.040742,5.254234)$del_num
plot_text_2_shape<-plot_text_2(LR_shape,9.2,a=7.040742,marker2_plot,siggene_shape_2$gene,N_del)


text(plot_text_2_shape$cm[1]+sub_rc[1]+20+20*(plot_text_2_shape$Pm[1]-1),
     plot_text_2_shape$LR[1]+sub_rc[4]+0.5,plot_text_2_shape$gene[1],adj=0,cex=0.6,family="Times New Roman1")

text(plot_text_2_shape$cm[2]+sub_rc[1]+20+20*(plot_text_2_shape$Pm[2]-1),
     plot_text_2_shape$LR[2]+sub_rc[4]+0.5,plot_text_2_shape$gene[2],adj=0,cex=0.6,family="Times New Roman1")


text(sub_rc[3]-70,sub_rc[2]-2,"C",cex=2,adj=0,family="Times New Roman")


#########################  D  petals  #################
sub_rc <-c(1020,25.5+20.5+20.5+20.5,1020+920, 5+20.5+20.5+20.5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

man_plot(LR_petals,68,marker2_plot,marker3_plot,a=12.227272,b=7.612247)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[3],sub_rc[4]+0.5+60*20/68/3*i,sub_rc[3]-8,sub_rc[4]+0.5+60*20/68/3*i,font=2)
  text(sub_rc[3]+40,sub_rc[4]+0.5+60*20/68/3*i,i*20,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_petals,12.227272),LR_petals,12.227272,7.612247)$del_num
plot_text_2_petals<-plot_text_2(LR_petals,68,a=12.227272,marker2_plot,siggene_petals_2$gene,N_del)

LR1<-plot_text_2_petals[which(plot_text_2_petals$Pm=="1"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:19){
  text_gene_plot2_100(LR1,i,20-0.8*i,0.6)
}


N_del<-gene_name1(signum_3(LR_petals,7.612247),LR_petals,12.227272,7.612247)$del_num
plot_text_3_petals<-plot_text_3(LR_petals,68,b=7.612247,marker3_plot,siggene_petals_3$gene,N_del)


LR1<-plot_text_3_petals[which(plot_text_3_petals$Pm=="1"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:9){
  text_gene_plot3_250(LR1,i,9-0.8*i,0.6)
}

text(sub_rc[3]-70,sub_rc[2]-2,"D",cex=2,adj=0,family="Times New Roman")


#############################  A  size ####################################
sub_rc <-c(100,25.5+20.5+20.5+20.5+20.5,1020, 5+20.5+20.5+20.5+20.5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

man_plot(LR_size,12.4,marker2_plot,marker3_plot,a=7.931283,b=7.334037)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+0.5+12*20/12.4/4*i,sub_rc[1]+8,sub_rc[4]+0.5+12*20/12.4/4*i,font=2)
  text(sub_rc[1]-40,sub_rc[4]+0.5+12*20/12.4/4*i,i*3,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_size,7.931283),LR_size,7.931283,7.334037)$del_num
plot_text_2_size<-plot_text_2(LR_size,12.4,a=7.931283,marker2_plot,siggene_size_2$gene,N_del)


LR1<-plot_text_2_size[which(plot_text_2_size$Pm=="2"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:4){
  text_gene_plot2_300(LR1,i,20-0.8*i,0.6)
}


LR1<-plot_text_2_size[which(plot_text_2_size$Pm=="6"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:4){
  text_gene_plot2_550(LR1,i,20-0.8*i,0.6)
}


text(plot_text_2_size$cm[9]+sub_rc[1]+20+20*(plot_text_2_size$Pm[9]-1),
     plot_text_2_size$LR[9]+sub_rc[4]+0.5,plot_text_2_size$gene[9],adj=0,cex=0.6,family="Times New Roman1")



N_del<-gene_name1(signum_3(LR_size,7.334037),LR_size,7.931283,7.334037)$del_num
plot_text_3_size<-plot_text_3(LR_size,12.4,b=7.334037,marker3_plot,siggene_size_3$gene,N_del)

text(plot_text_3_size$cm[1]+sub_rc[1]+20+20*(plot_text_3_size$Pm[1]-1),
     plot_text_3_size$LR[1]+sub_rc[4]+0.5,plot_text_3_size$gene[1],adj=0,cex=0.6,family="Times New Roman1")


text(plot_text_3_size$cm[2]+sub_rc[1]+20+20*(plot_text_3_size$Pm[2]-1),
     plot_text_3_size$LR[2]+sub_rc[4]+0.5,plot_text_3_size$gene[2],adj=0,cex=0.6,family="Times New Roman1")

text(plot_text_3_size$cm[3]+sub_rc[1]+20+20*(plot_text_3_size$Pm[3]-1),
     plot_text_3_size$LR[3]+sub_rc[4]+0.5,plot_text_3_size$gene[3],adj=0,cex=0.6,family="Times New Roman1")

text(sub_rc[1]+10,sub_rc[2]-2,"A",cex=2,adj=0,family="Times New Roman")


#########################  B  pedicels  #################
sub_rc <-c(1020,25.5+20.5+20.5+20.5+20.5,1020+920, 5+20.5+20.5+20.5+20.5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

man_plot(LR_pedicels,10,marker2_plot,marker3_plot,a=6.575807,b=6.402251)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[3],sub_rc[4]+0.5+9*20/10/3*i,sub_rc[3]-8,sub_rc[4]+0.5+9*20/10/3*i,font=2)
  text(sub_rc[3]+40,sub_rc[4]+0.5+9*20/10/3*i,i*3,cex=1.2,font=1,family="Times New Roman")
}


N_del<-gene_name1(signum_2(LR_pedicels,6.575807),LR_pedicels,6.575807,6.402251)$del_num
plot_text_2_pedicels<-plot_text_2(LR_pedicels,10,a=6.575807,marker2_plot,siggene_pedicels_2$gene,N_del)


text(plot_text_2_pedicels$cm[1]+sub_rc[1]+20+20*(plot_text_2_pedicels$Pm[1]-1),
     plot_text_2_pedicels$LR[1]+sub_rc[4]+0.5,plot_text_2_pedicels$gene[1],adj=1,cex=0.6,family="Times New Roman1")


N_del<-gene_name1(signum_3(LR_pedicels,6.402251),LR_pedicels,6.575807,6.402251)$del_num
plot_text_3_pedicels<-plot_text_3(LR_pedicels,10,b=6.402251,marker3_plot,siggene_pedicels_3$gene,N_del)

text(plot_text_3_pedicels$cm[1]+sub_rc[1]+10+20*(plot_text_3_pedicels$Pm[1]-1),
     plot_text_3_pedicels$LR[1]+sub_rc[4]+1.5,plot_text_3_pedicels$gene[1],adj=1,cex=0.6,family="Times New Roman1")


segments(plot_text_3_pedicels$cm[1]+sub_rc[1]+20+20*(plot_text_3_pedicels$Pm[1]-1),
         plot_text_3_pedicels$LR[1]+sub_rc[4]+0.5,
         plot_text_3_pedicels$cm[1]+sub_rc[1]+10+20*(plot_text_3_pedicels$Pm[1]-1),
         plot_text_3_pedicels$LR[1]+sub_rc[4]+1.5,lwd=0.5)
text(sub_rc[1]+10,sub_rc[2]-2,"B",cex=2,adj=0,family="Times New Roman")


dev.off()
