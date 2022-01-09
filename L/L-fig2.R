
pdf("2020-1/L-fig2.pdf",width=10,height =3)


par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,4300),ylim=c(0, 28.5))

#############################  A  ####################################
sub_rc <-c(200,25.5,2280, 5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

man_plot(LR_shape,12,marker2_plot,marker3_plot,a=7.101772,b=6.717516)


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+0.5+5*i,sub_rc[1]+15,sub_rc[4]+0.5+5*i,font=2)
  text(sub_rc[1]-60,sub_rc[4]+0.5+5*i,i*3,cex=1.2,font=1,family="Times New Roman")
}


A<-c(1:8)
for(i in 1:8){
  #segments(sub_rc[1]+50+bpMidVec[i]+20*(i-1),sub_rc[4],sub_rc[1]+50+bpMidVec[i]+20*(i-1),sub_rc[4]+0.4,font=2)
  text(sub_rc[1]+50+bpMidVec[i]+20*(i-1),sub_rc[4]-1,A[i],cex=1.2,font=1,family="Times New Roman")
}

#segments(sub_rc[1]+50+bpMidVec[9]+20*8,sub_rc[4],sub_rc[1]+50+bpMidVec[9]+20*8,sub_rc[4]+0.4,font=2)
text(sub_rc[1]+50+bpMidVec[9]+20*8,sub_rc[4]-1,"scaffold",cex=0.8,font=1,adj=0,family="Times New Roman")



N_del<-gene_name1(signum_2(LR_shape,7.101772),LR_shape,7.101772,6.717516)$del_num
plot_text_2_shape<-plot_text_2(LR_shape,12,a=7.101772,marker2_plot,siggene_shape_2$gene,N_del)

LR1<-plot_text_2_shape[which(plot_text_2_shape$Pm=="1"),]
LR1<-LR1[which(duplicated(LR1$gene, fromLast=TRUE)==FALSE),]
rownames(LR1)<-NULL
LR1<-LR1[order(LR1$LR, decreasing = T),]

for(i in 1:15){
  text_gene_plot2_300(LR1,i,20.5-1*i,0.8)
}


text(plot_text_2_shape$cm[19]+sub_rc[1]+50+20*(plot_text_2_shape$Pm[19]-1),
     plot_text_2_shape$LR[19]+sub_rc[4]+0.5,plot_text_2_shape$gene[19],adj=0,cex=0.8,family="Times New Roman1")

text(plot_text_2_shape$cm[20]+sub_rc[1]+50+20*(plot_text_2_shape$Pm[20]-1),
     plot_text_2_shape$LR[20]+sub_rc[4]+0.5,plot_text_2_shape$gene[20],adj=0,cex=0.8,family="Times New Roman1")


N_del<-gene_name1(signum_3(LR_shape,6.717516),LR_shape,7.101772,6.717516)$del_num
plot_text_3_shape<-plot_text_3(LR_shape,12,b=6.717516,marker3_plot,siggene_shape_3$gene,N_del)

text(plot_text_3_shape$cm[1]+sub_rc[1]+650+20*(plot_text_3_shape$Pm[1]-1),
     plot_text_3_shape$LR[1]+sub_rc[4]+0.5,plot_text_3_shape$gene[1],adj=0,cex=0.8,family="Times New Roman1")

segments(plot_text_3_shape$cm[1]+sub_rc[1]+50+20*(plot_text_3_shape$Pm[1]-1),
         plot_text_3_shape$LR[1]+sub_rc[4]+0.5,
         plot_text_3_shape$cm[1]+sub_rc[1]+650+20*(plot_text_3_shape$Pm[1]-1),
         plot_text_3_shape$LR[1]+sub_rc[4]+0.5,lwd=0.5) 

text(sub_rc[1]-220,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression(-log(p-value)),srt=90,cex=1.6,family="Times New Roman")
text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-4.5,"Chromosome",cex=2,family="Times New Roman")
text(sub_rc[1],sub_rc[2]+3,"A",cex=2,adj=0,family="Times New Roman")


########################################

sub_rc <-c(200+2300,25.5,2280+1020, 5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


gene="Pm001721"

marker_num<-siggene_shape_2$Marker_ID[which(siggene_shape_2$gene==gene)]


n1<-which(L_marker_st$marker[,3][-1]==marker_num)

geno_marker(dat_L1_st$geno[,n1])
marker_pheno_shape_1<-pheno_geno(dat_L1_st$geno[,n1],dat_L1_st$pheno[,4])

ratio1<-length(which(marker_pheno_shape_1[[1]]==1))/length(marker_pheno_shape_1[[1]])

ratio2<-length(which(marker_pheno_shape_1[[2]]==1))/length(marker_pheno_shape_1[[2]])




rect(sub_rc[1]+140,sub_rc[2],sub_rc[1]+140+190,sub_rc[2]-ratio1*20.5,border="#E41A1C",col="#E41A1C",lwd=1)

rect(sub_rc[1]+140,sub_rc[2]-ratio1*20.5,sub_rc[1]+140+190,sub_rc[4],border="#FFFF33",col="#FFFF33",lwd=1)

rect(sub_rc[1]+140+140+190,sub_rc[2],sub_rc[1]+140+190+140+190,sub_rc[2]-ratio2*20.5,border="#E41A1C",col="#E41A1C",lwd=1)

rect(sub_rc[1]+140+140+190,sub_rc[2]-ratio2*20.5,sub_rc[1]+140+190+140+190,sub_rc[4],border="#FFFF33",col="#FFFF33",lwd=1)

text(sub_rc[1]+140+190/2,sub_rc[4]-1.2,"ll",cex=1.5,adj=0.5,family="Times New Roman1")
text(sub_rc[1]+140+190/2+140+190,sub_rc[4]-1.2,"lm",cex=1.5,adj=0.5,family="Times New Roman1")

text((sub_rc[3]-sub_rc[1])/2+sub_rc[1],sub_rc[2]+1.2,"Pm001721",cex=1.2,adj=0.5,family="Times New Roman1")

text(sub_rc[1],sub_rc[2]+3,"B",cex=2,adj=0,family="Times New Roman")

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)




sub_rc <-c(200+2300+800,25.5,2280+1020+800, 5)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


gene="Pm000784"

marker_num<-siggene_shape_3$Marker_ID[which(siggene_shape_3$gene==gene)]


n1<-which(L_marker_st$marker[,3][-1]==marker_num)

geno_marker(dat_L1_st$geno[,n1])
marker_pheno_shape_1<-pheno_geno(dat_L1_st$geno[,n1],dat_L1_st$pheno[,4])

ratio1<-length(which(marker_pheno_shape_1[[1]]==1))/length(marker_pheno_shape_1[[1]])

ratio2<-length(which(marker_pheno_shape_1[[2]]==1))/length(marker_pheno_shape_1[[2]])

ratio3<-length(which(marker_pheno_shape_1[[3]]==1))/length(marker_pheno_shape_1[[3]])


rect(sub_rc[1]+80,sub_rc[2],sub_rc[1]+80+160,sub_rc[2]-ratio1*20.5,border="#E41A1C",col="#E41A1C",lwd=1)

rect(sub_rc[1]+80,sub_rc[2]-ratio1*20.5,sub_rc[1]+80+160,sub_rc[4],border="#FFFF33",col="#FFFF33",lwd=1)

rect(sub_rc[1]+80+80+160,sub_rc[2],sub_rc[1]+80+160+80+160,sub_rc[2]-ratio2*20.5,border="#E41A1C",col="#E41A1C",lwd=1)

rect(sub_rc[1]+80+80+160,sub_rc[2]-ratio2*20.5,sub_rc[1]+80+160+80+160,sub_rc[4],border="#FFFF33",col="#FFFF33",lwd=1)

rect(sub_rc[1]+80+80+160+80+160,sub_rc[2],sub_rc[1]+80+160+80+160+80+160,sub_rc[2]-ratio3*20.5,border="#E41A1C",col="#E41A1C",lwd=1)

rect(sub_rc[1]+80+80+160+80+160,sub_rc[2]-ratio3*20.5,sub_rc[1]+80+160+80+160+80+160,sub_rc[4],border="#FFFF33",col="#FFFF33",lwd=1)


text(sub_rc[1]+80+160/2,sub_rc[4]-1.2,"hh",cex=1.5,adj=0.5,family="Times New Roman1")
text(sub_rc[1]+80+160/2+80+160,sub_rc[4]-1.2,"hk",cex=1.5,adj=0.5,family="Times New Roman1")
text(sub_rc[1]+80+160/2+80+160+80+160,sub_rc[4]-1.2,"kk",cex=1.5,adj=0.5,family="Times New Roman1")

text((sub_rc[3]-sub_rc[1])/2+sub_rc[1],sub_rc[2]+1.2,"Pm000784",cex=1.2,adj=0.5,family="Times New Roman1")

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[3],sub_rc[4]+20.5/5*i,sub_rc[3]-15,sub_rc[4]++20.5/5*i,font=2)
  text(sub_rc[3]+80,sub_rc[4]++20.5/5*i,20*i,cex=1.2,font=1,family="Times New Roman")
}

text(sub_rc[3]+220,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],"Percentage (%)",srt=90,cex=1.2,family="Times New Roman")


text(sub_rc[1],sub_rc[4]-4.5,"Genotype",cex=2,adj=0.5,family="Times New Roman")

dev.off()



