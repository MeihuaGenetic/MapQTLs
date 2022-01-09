
pdf("2020-1/Y.pdf",width=10,height =14.4)


par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(40,2700),ylim=c(4, 28.5+20.5+20.5+20.5+20.5+20.5-5))


#############################  K   daimeter  ####################################

sub_rc <-c(180,25.5,1380, 5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


man_plot(LR_daimeter,8,marker2_plot,marker3_plot,6.568583,5.546653)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[1],sub_rc[4]+0.5+8*20/8/4*i,sub_rc[1]+8,sub_rc[4]+0.5+8*20/8/4*i,font=2)
  text(sub_rc[1]-40,sub_rc[4]+0.5+8*20/8/4*i,i*2,cex=1.2,font=1,family="Times New Roman")
}

A<-c(1:8)
for(i in 1:8){
  #segments(sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4],sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4]+0.4,font=2)
  text(sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4]-1,A[i],cex=1.2,font=1,family="Times New Roman")
}

#segments(sub_rc[1]+20+bpMidVec[9]+20*8,sub_rc[4],sub_rc[1]+20+bpMidVec[9]+20*8,sub_rc[4]+0.4,font=2)
text(sub_rc[1]+20+bpMidVec[9]+20*8,sub_rc[4]-1,"scaffold",cex=0.8,font=1,adj=0,family="Times New Roman")

N_del<-gene_name1(signum_2(LR_daimeter,6.568583),LR_daimeter,6.568583,5.546653)$del_num
plot_text_2_daimeter<-plot_text_2(LR_daimeter,8,a=6.568583,marker2_plot,siggene_daimeter_2$gene,N_del)

N_del<-gene_name1(signum_3(LR_daimeter,5.546653),LR_daimeter,6.568583,5.546653)$del_num
plot_text_3_daimeter<-plot_text_3(LR_daimeter,8,b=5.546653,marker3_plot,siggene_daimeter_3$gene,N_del)

text(plot_text_3_daimeter$cm[1]+sub_rc[1]+20+20*(plot_text_3_daimeter$Pm[1]-1),
     plot_text_3_daimeter$LR[1]+sub_rc[4]+0.5,plot_text_3_daimeter$gene[1],adj=0,cex=0.8,family="Times New Roman1")

text(plot_text_3_daimeter$cm[2]+sub_rc[1]+20+20*(plot_text_3_daimeter$Pm[2]-1),
     plot_text_3_daimeter$LR[2]+sub_rc[4]+0.5,plot_text_3_daimeter$gene[2],adj=0,cex=0.8,family="Times New Roman1")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-3.5,"Chromosome",cex=2,family="Times New Roman")
text(sub_rc[1]+10,sub_rc[2]-2,"K",cex=2,adj=0,family="Times New Roman")



############################# I   period ####################################

sub_rc <-c(180,25.5+20.5,1380, 5+20.5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


man_plot(LR_period,9,marker2_plot,marker3_plot,7.227870,7.177984)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+0.5+8*20/9/4*i,sub_rc[1]+8,sub_rc[4]+0.5+8*20/9/4*i,font=2)
  text(sub_rc[1]-40,sub_rc[4]+0.5+8*20/9/4*i,i*2,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_period,7.227870),LR_period,7.227870,7.177984)$del_num
plot_text_2_period<-plot_text_2(LR_period,9,a=7.227870,marker2_plot,siggene_period_2$gene,N_del)

text(plot_text_2_period$cm[1]+sub_rc[1]+20+20*(plot_text_2_period$Pm[1]-1),
     plot_text_2_period$LR[1]+sub_rc[4]+0.5,plot_text_2_period$gene[1],adj=0,cex=0.8,family="Times New Roman1")

text(plot_text_2_period$cm[2]+sub_rc[1]+20+20*(plot_text_2_period$Pm[2]-1),
     plot_text_2_period$LR[2]+sub_rc[4]+0.5,plot_text_2_period$gene[2],adj=0,cex=0.8,family="Times New Roman1")

text(plot_text_2_period$cm[3]+sub_rc[1]+20+20*(plot_text_2_period$Pm[3]-1),
     plot_text_2_period$LR[3]+sub_rc[4]+0.5,plot_text_2_period$gene[3],adj=0,cex=0.8,family="Times New Roman1")

text(sub_rc[1]+10,sub_rc[2]-2,"I",cex=2,adj=0,family="Times New Roman")


############################# J   length ####################################

sub_rc <-c(1380,25.5+20.5,1380+1200,5+20.5)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


man_plot(LR_length,7,marker2_plot,marker3_plot,a=6.510982,b=6.000698)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:3){
  segments(sub_rc[3],sub_rc[4]+0.5+6*20/7/3*i,sub_rc[3]-8,sub_rc[4]+0.5+6*20/7/3*i,font=2)
  text(sub_rc[3]+40,sub_rc[4]+0.5+6*20/7/3*i,i*2,cex=1.2,font=1,family="Times New Roman")
}

text(sub_rc[1]+10,sub_rc[2]-2,"J",cex=2,adj=0,family="Times New Roman")

A<-c(1:8)
for(i in 1:8){
  #segments(sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4],sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4]+0.4,font=2)
  text(sub_rc[1]+20+bpMidVec[i]+20*(i-1),sub_rc[4]-1,A[i],cex=1.2,font=1,family="Times New Roman")
}

#segments(sub_rc[1]+20+bpMidVec[9]+20*8,sub_rc[4],sub_rc[1]+20+bpMidVec[9]+20*8,sub_rc[4]+0.4,font=2)
text(sub_rc[1]+20+bpMidVec[9]+20*8,sub_rc[4]-1,"scaffold",cex=0.8,font=1,adj=0,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-3.5,"Chromosome",cex=2,family="Times New Roman")


############################# G   center####################################

sub_rc <-c(180,25.5+20.5+20.5,1380, 5+20.5+20.5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
man_plot(LR_center,8,marker2_plot,marker3_plot,5.771643,5.997276)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[1],sub_rc[4]+0.5+8*20/8/4*i,sub_rc[1]+8,sub_rc[4]+0.5+8*20/8/4*i,font=2)
  text(sub_rc[1]-40,sub_rc[4]+0.5+8*20/8/4*i,i*2,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_center,5.771643),LR_center,5.771643,5.997276)$del_num
plot_text_2_center<-plot_text_2(LR_center,8,a=5.771643,marker2_plot,siggene_center_2$gene,N_del)

text(plot_text_2_center$cm[1]+sub_rc[1]+20+20*(plot_text_2_center$Pm[1]-1),
     plot_text_2_center$LR[1]+sub_rc[4]+0.5,plot_text_2_center$gene[1],adj=0,cex=0.8,family="Times New Roman1")

text(plot_text_2_center$cm[2]+sub_rc[1]+20+20*(plot_text_2_center$Pm[2]-1),
     plot_text_2_center$LR[2]+sub_rc[4]+0.5,plot_text_2_center$gene[2],adj=1,cex=0.8,family="Times New Roman1")

text(plot_text_2_center$cm[3]+sub_rc[1]+20+20*(plot_text_2_center$Pm[3]-1),
     plot_text_2_center$LR[3]+sub_rc[4]+0.5,plot_text_2_center$gene[3],adj=1,cex=0.8,family="Times New Roman1")

text(sub_rc[1]+10,sub_rc[2]-2,"G",cex=2,adj=0,family="Times New Roman")


############################# H   colour ####################################

sub_rc <-c(1380,25.5+20.5+20.5,1380+1200,5+20.5+20.5)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
man_plot(LR_colour,10,marker2_plot,marker3_plot,5.897320,6.181876)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[3],sub_rc[4]+0.5+9*20/10/3*i,sub_rc[3]-8,sub_rc[4]+0.5+9*20/10/3*i,font=2)
  text(sub_rc[3]+40,sub_rc[4]+0.5+9*20/10/3*i,i*3,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_colour,5.897320),LR_colour,5.897320,6.181876)$del_num
plot_text_2_colour<-plot_text_2(LR_colour,10,a=5.897320,marker2_plot,siggene_colour_2$gene,N_del)

text(plot_text_2_colour$cm[1]+sub_rc[1]+20+20*(plot_text_2_colour$Pm[1]-1),
     plot_text_2_colour$LR[1]+sub_rc[4]+0.5,plot_text_2_colour$gene[1],adj=0,cex=0.8,family="Times New Roman1")

text(plot_text_2_colour$cm[2]+sub_rc[1]+20+20*(plot_text_2_colour$Pm[2]-1),
     plot_text_2_colour$LR[2]+sub_rc[4]+0.5,plot_text_2_colour$gene[2],adj=0,cex=0.8,family="Times New Roman1")

text(plot_text_2_colour$cm[3]+sub_rc[1]+20+20*(plot_text_2_colour$Pm[3]-1),
     plot_text_2_colour$LR[3]+sub_rc[4]+0.5,plot_text_2_colour$gene[3],adj=0,cex=0.8,family="Times New Roman1")

text(sub_rc[1]+10,sub_rc[2]-2,"H",cex=2,adj=0,family="Times New Roman")


############################# E   center####################################

sub_rc <-c(180,25.5+20.5+20.5+20.5,1380, 5+20.5+20.5+20.5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
man_plot(LR_pistils,7.4,marker2_plot,marker3_plot,6.557408 ,5.101051)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[1],sub_rc[4]+0.5+6*20/7.4/3*i,sub_rc[1]+8,sub_rc[4]+0.5+6*20/7.4/3*i,font=2)
  text(sub_rc[1]-40,sub_rc[4]+0.5+6*20/7.4/3*i,i*2,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_3(LR_pistils,5.101051),LR_bud,6.557408 ,5.101051)$del_num
plot_text_3_pistils<-plot_text_3(LR_pistils,7.4,b=5.101051,marker3_plot,siggene_pistils_3$gene,N_del)

text(plot_text_3_pistils$cm[1]+sub_rc[1]+20+20*(plot_text_3_pistils$Pm[1]-1),
     plot_text_3_pistils$LR[1]+sub_rc[4]+0.5,plot_text_3_pistils$gene[1],adj=0,cex=0.8,family="Times New Roman1")

text(sub_rc[1]+10,sub_rc[2]-2,"E",cex=2,adj=0,family="Times New Roman")

text(sub_rc[1]-140,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression(-log(p-value)),srt=90,cex=1.6,family="Times New Roman")

############################# F   colour ####################################

sub_rc <-c(1380,25.5+20.5+20.5+20.5,1380+1200,5+20.5+20.5+20.5)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
man_plot(LR_bud,6,marker2_plot,marker3_plot,4.757960, 4.672774)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:2){
  segments(sub_rc[3],sub_rc[4]+0.5+6*20/6/3*i,sub_rc[3]-8,sub_rc[4]+0.5+6*20/6/3*i,font=2)
  text(sub_rc[3]+40,sub_rc[4]+0.5+6*20/6/3*i,i*2,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_bud,4.757960),LR_bud,4.757960, 4.672774)$del_num
plot_text_2_bud<-plot_text_2(LR_bud,6,a=4.757960,marker2_plot,siggene_bud_2$gene,N_del)

text(plot_text_2_bud$cm[1]+sub_rc[1]+20+20*(plot_text_2_bud$Pm[1]-1),
     plot_text_2_bud$LR[1]+sub_rc[4]+0.5,plot_text_2_bud$gene[1],adj=0,cex=0.8,family="Times New Roman1")

text(plot_text_2_bud$cm[2]+sub_rc[1]+20+20*(plot_text_2_bud$Pm[2]-1),
       plot_text_2_bud$LR[2]+sub_rc[4]+0.5,plot_text_2_bud$gene[2],adj=1,cex=0.8,family="Times New Roman1")


text(sub_rc[1]+10,sub_rc[2]-2,"F",cex=2,adj=0,family="Times New Roman")

text(sub_rc[3]+140,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression(-log(p-value)),srt=90,cex=1.6,family="Times New Roman")

############################# C   center####################################

sub_rc <-c(180,25.5+20.5+20.5+20.5+20.5,1380, 5+20.5+20.5+20.5+20.5)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
man_plot(LR_shape,11,marker2_plot,marker3_plot,6.563671,6.136940)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[1],sub_rc[4]+0.5+9*20/11/3*i,sub_rc[1]+8,sub_rc[4]+0.5+9*20/11/3*i,font=2)
  text(sub_rc[1]-40,sub_rc[4]+0.5+9*20/11/3*i,i*3,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_shape,6.563671),LR_shape,6.563671,6.136940)$del_num
plot_text_2_shape<-plot_text_2(LR_shape,11,a=6.563671,marker2_plot,siggene_shape_2$gene,N_del)


LR1<-plot_text_2_shape[which(plot_text_2_shape$Pm=="6"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:6){
  text_gene_plot2_750(LR1,i,20-1.5*i,0.8)
}


text(sub_rc[1]+10,sub_rc[2]-2,"C",cex=2,adj=0,family="Times New Roman")


############################# D   petals ####################################

sub_rc <-c(1380,25.5+20.5+20.5+20.5+20.5,1380+1200,5+20.5+20.5+20.5+20.5)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
man_plot(LR_petals,7.4,marker2_plot,marker3_plot,7.029034 ,4.920261)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:3){
  segments(sub_rc[3],sub_rc[4]+0.5+6*20/7.4/3*i,sub_rc[3]-8,sub_rc[4]+0.5+6*20/7.4/3*i,font=2)
  text(sub_rc[3]+40,sub_rc[4]+0.5+6*20/7.4/3*i,i*2,cex=1.2,font=1,family="Times New Roman")
}


text(sub_rc[1]+10,sub_rc[2]-2,"D",cex=2,adj=0,family="Times New Roman")

############################# A   size####################################

sub_rc <-c(180,25.5+20.5+20.5+20.5+20.5+20.5,1380, 5+20.5+20.5+20.5+20.5+20.5)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
man_plot(LR_size,8,marker2_plot,marker3_plot,6.827272 ,6.072412)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+0.5+8*20/8/4*i,sub_rc[1]+8,sub_rc[4]+0.5+8*20/8/4*i,font=2)
  text(sub_rc[1]-40,sub_rc[4]+0.5+8*20/8/4*i,i*2,cex=1.2,font=1,family="Times New Roman")
}


N_del<-gene_name1(signum_3(LR_size,6.072412),LR_size,6.827272 ,6.072412)$del_num
plot_text_3_size<-plot_text_3(LR_size,8,b=6.072412,marker3_plot,siggene_size_3$gene,N_del)

text(plot_text_3_size$cm[1]+sub_rc[1]+20+20*(plot_text_3_size$Pm[1]-1),
     plot_text_3_size$LR[1]+sub_rc[4]+0.5,plot_text_3_size$gene[1],adj=0,cex=0.8,family="Times New Roman1")

text(plot_text_3_size$cm[2]+sub_rc[1]+20+20*(plot_text_3_size$Pm[2]-1),
     plot_text_3_size$LR[2]+sub_rc[4]+0.5,plot_text_3_size$gene[2],adj=0,cex=0.8,family="Times New Roman1")

text(sub_rc[1]+10,sub_rc[2]-2,"A",cex=2,adj=0,family="Times New Roman")

############################# B   size####################################

sub_rc <-c(1380,25.5+20.5+20.5+20.5+20.5+20.5,1380+1200,5+20.5+20.5+20.5+20.5+20.5)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
man_plot(LR_pedicels,8,marker2_plot,marker3_plot,6.618790,6.041480)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:4){
  segments(sub_rc[3],sub_rc[4]+0.5+8*20/8/4*i,sub_rc[3]-8,sub_rc[4]+0.5+8*20/8/4*i,font=2)
  text(sub_rc[3]+40,sub_rc[4]+0.5+8*20/8/4*i,i*2,cex=1.2,font=1,family="Times New Roman")
}

N_del<-gene_name1(signum_2(LR_pedicels,6.618790),LR_pedicels,6.618790,6.041480)$del_num
plot_text_2_pedicels<-plot_text_2(LR_pedicels,8,a=6.618790,marker2_plot,siggene_pedicels_2$gene,N_del)

#plot_text_3_pedicels<-plot_text_3(LR_pedicels,b=7.726914,marker3_plot,siggene_pedicels_3$gene)

text(plot_text_2_pedicels$cm[1]+sub_rc[1]+20+20*(plot_text_2_pedicels$Pm[1]-1),
     plot_text_2_pedicels$LR[1]+sub_rc[4]+0.5,plot_text_2_pedicels$gene[1],adj=0,cex=0.8,family="Times New Roman1")


text(sub_rc[1]+10,sub_rc[2]-2,"B",cex=2,adj=0,family="Times New Roman")

dev.off()

