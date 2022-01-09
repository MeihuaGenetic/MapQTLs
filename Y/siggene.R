
signum_anno_1<-function(LR_calyx,a,b){
  marker_table <-  Y_marker_st$marker[-1,1:2]
  
  info_1<-cbind(marker_table,LR_calyx)
  
  marker_2<-info_1[geno_type_st$marker2,]
  marker_3<-info_1[geno_type_st$marker3,]
  
  num2<-geno_type_st$marker2[which(marker_2[,3]>a)]
  num3<-geno_type_st$marker3[which(marker_3[,3]>b)]
  num<-sort(c(num2,num3))
  num
}



marker_name_anno_1<-function(LR_calyx,a,b){
  marker_table2 <- Y_marker_st$marker[-1,1:2]
  signum<-signum_anno_1(LR_calyx,a,b)
  marker_table2 <-  marker_table2[signum,]
  na_num<-which(marker_table2[,2]=="None")
  if(length(na_num)==0){
    marker_table1 <-marker_table2
  }else{
    marker_table1 <-marker_table2[-na_num,]
  }
  
  name1<-as.matrix(gsub("Pm","Pa",marker_table1[,1]),ncol=1)
  marker_table<-cbind(name1,as.numeric(as.character(marker_table1[,2])))
  
  plum_gff<-read.table("plum.genome.gff")
  plum_gff1<-plum_gff[which(plum_gff[,3]=="mRNA"),]
  
  marker_name<-names(table(marker_table[,1]))
  plum_name<-names(table(plum_gff1[,1]))
  
  marker<-c()
  for(i in 1:dim(marker_table)[1]){
    A<-marker_table[i,1]
    if(length(which(A==plum_name))==1){
      marker<-rbind(marker,marker_table[i,])
    }
  }
  marker<-as.data.frame(marker)
  colnames(marker)<-c("PmID","PmPosition")
  marker
}

a=7.260681
b=6.915640


gene_name_anno_1<-function(LR_calyx,a,b){
  marker_name<-marker_name_anno_1(LR_calyx,a,b)
  plum_gff<-read.table("plum.genome.gff")
  plum_gff1<-plum_gff[which(plum_gff[,3]=="mRNA"),]  
  A5<-c()
  for(i in 1:dim(marker_name)[1]){
    Pm<-as.character(marker_name$PmID[i])
    AA<-plum_gff1[which(Pm==plum_gff1[,1]),]
    
    if(dim(AA)[1]!=0){
      pos<-as.numeric(as.character(marker_name$PmPosition[i]))
      
      gene_1<-c()
      if(max(as.numeric(AA[,5]))<pos){
        gene_1<-"None"
      }else{
        if(min(as.numeric(AA[,4]))>pos){
          gene_1<-"None"
        }else{
          num1<-max(which(as.numeric(AA[,4])<pos))
          num2<-min(which(as.numeric(AA[,5])>pos))
          if(num1!=num2){
            gene_1<-"None"
          }else{
            num<-num1
            gene_1<-substr(AA[num,9],4,11)
          }
        }
      }
      
      A5<-rbind(A5,gene_1)
    }
  }
  A6<-cbind(marker_name,A5)
  A6<-as.data.frame(A6)
  colnames(A6)<-c("PmID","PmPosition","gene")
  rownames(A6)<-NULL
  A6
}



KEGG_anno1<-function(gene_name){
  KEGG_anno<-matrix(NA,dim(gene_name)[1],1) 
  for(i in 1:dim(gene_name)[1]){
    pm<-as.character(gene_name[i,3])
    num<-which(KEGG[,1]==pm)
    if(length(num)!=0){
      KEGG_anno[i,]<-as.character(KEGG[num,16])
    }
  }
  KEGG_anno
}

GO_anno1<-function(gene_name){
  GO_anno1<-matrix(NA,dim(gene_name)[1],1) 
  for(i in 1:dim(gene_name)[1]){
    pm<-as.character(gene_name[i,3])
    num<-which(GO[,1]==pm)
    if(length(num)!=0){
      GO1<-GO[num,3:12]
      GO1<-na.omit(c(as.matrix(GO1,nrow=1)))
      num2<-length(GO1)
      if(num2==1)
        GO2<-GO1
      if(num2==2)
        GO2<-paste(GO1[1],GO1[2], sep = " ")
      if(num2==3)
        GO2<-paste(GO1[1],GO1[2],GO1[3], sep = " ")
      if(num2==4)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4], sep = " ")
      if(num2==5)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5], sep = " ")
      if(num2==6)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6], sep = " ")
      if(num2==7)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7], sep = " ")
      if(num2==8)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8], sep = " ")
      if(num2==9)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9], sep = " ")
      if(num2==10)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], sep = " ")
      if(num2==11)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],sep = " ")
      if(num2==12)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],sep = " ")
      if(num2==13)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],GO1[13],sep = " ")
      if(num2==14)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],GO1[13],GO1[14],sep = " ")
      if(num2==15)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],GO1[13],GO1[14],GO1[15],sep = " ")
      if(num2==16)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],GO1[13],GO1[14],GO1[15],GO1[16],sep = " ")
      if(num2==17)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],GO1[13],GO1[14],GO1[15],GO1[16],GO1[17],sep = " ")
      GO_anno1[i,]<-GO2
    }
  }
  GO_anno2<-matrix(NA,dim(gene_name)[1],1) 
  for(i in 1:dim(gene_name)[1]){
    pm<-as.character(gene_name[i,3])
    num<-which(GO[,1]==pm)
    if(length(num)!=0){
      GO_anno2[i,]<-GO[num,2]
    }
  }
  GO_anno<-cbind(GO_anno2,GO_anno1)
  GO_anno
}



GENE_anno1<-function(gene_name){
  GO_anno1<-matrix(NA,dim(gene_name)[1],1) 
  for(i in 1:dim(gene_name)[1]){
    pm<-as.character(gene_name[i,3])
    num<-which(GENE[,1]==pm)
    if(length(num)!=0){
      GO1<-GENE[num,3:12]
      GO1<-na.omit(c(as.matrix(GO1,nrow=1)))
      num2<-length(GO1)
      if(num2==1)
        GO2<-GO1
      if(num2==2)
        GO2<-paste(GO1[1],GO1[2], sep = " ")
      if(num2==3)
        GO2<-paste(GO1[1],GO1[2],GO1[3], sep = " ")
      if(num2==4)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4], sep = " ")
      if(num2==5)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5], sep = " ")
      if(num2==6)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6], sep = " ")
      if(num2==7)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7], sep = " ")
      if(num2==8)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8], sep = " ")
      if(num2==9)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9], sep = " ")
      if(num2==10)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], sep = " ")
      if(num2==11)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],sep = " ")
      if(num2==12)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],sep = " ")
      if(num2==13)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],GO1[13],sep = " ")
      if(num2==14)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],GO1[13],GO1[14],sep = " ")
      if(num2==15)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],GO1[13],GO1[14],GO1[15],sep = " ")
      if(num2==16)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],GO1[13],GO1[14],GO1[15],GO1[16],sep = " ")
      if(num2==17)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],GO1[13],GO1[14],GO1[15],GO1[16],GO1[17],sep = " ")
      if(num2==17)
        GO2<-paste(GO1[1],GO1[2],GO1[3],GO1[4],GO1[5],GO1[6],GO1[7],GO1[8],GO1[9],GO1[10], 
                   GO1[11],GO1[12],GO1[13],GO1[14],GO1[15],GO1[16],GO1[17],GO1[18],sep = " ")
      GO_anno1[i,]<-GO2
    }
  }
  GO_anno2<-matrix(NA,dim(gene_name)[1],1) 
  for(i in 1:dim(gene_name)[1]){
    pm<-as.character(gene_name[i,3])
    num<-which(GENE[,1]==pm)
    if(length(num)!=0){
      GO_anno2[i,]<-GENE[num,2]
    }
  }
  GO_anno<-cbind(GO_anno2,GO_anno1)
  GO_anno
}



gene_anno_colour<-gene_name_anno_1(LR_colour,5.897320,6.181876)
gene_anno_bud<-gene_name_anno_1(LR_bud,4.757960, 4.672774)
gene_anno_shape<-gene_name_anno_1(LR_shape,6.563671,6.136940)
gene_anno_period<-gene_name_anno_1(LR_period,7.227870,7.177984)
gene_anno_center<-gene_name_anno_1(LR_center,5.771643,5.997276)
gene_anno_pistils<-gene_name_anno_1(LR_pistils,6.557408 ,5.101051)
gene_anno_size<-gene_name_anno_1(LR_size,6.827272 ,6.072412)
gene_anno_petals<-gene_name_anno_1(LR_petals,7.029034 ,4.920261)
gene_anno_pedicels<-gene_name_anno_1(LR_pedicels,6.618790,6.041480)
gene_anno_length<-gene_name_anno_1(LR_length,6.510982,6.000698)
gene_anno_daimeter<-gene_name_anno_1(LR_daimeter,6.568583,5.546653)


write.csv(gene_anno_colour,file="2020-1/siggene/Y-colour.csv",row.names = FALSE)
write.csv(gene_anno_bud,file="2020-1/siggene/Y-bud.csv",row.names = FALSE)
write.csv(gene_anno_shape,file="2020-1/siggene/Y-shape.csv",row.names = FALSE)
write.csv(gene_anno_period,file="2020-1/siggene/Y-period.csv",row.names = FALSE)
write.csv(gene_anno_pistils,file="2020-1/siggene/Y-pistils.csv",row.names = FALSE)
write.csv(gene_anno_size,file="2020-1/siggene/Y-size.csv",row.names = FALSE)
write.csv(gene_anno_petals,file="2020-1/siggene/Y-petals.csv",row.names = FALSE)
write.csv(gene_anno_pedicels,file="2020-1/siggene/Y-pedicels.csv",row.names = FALSE)
write.csv(gene_anno_length,file="2020-1/siggene/Y-length.csv",row.names = FALSE)
write.csv(gene_anno_daimeter,file="2020-1/siggene/Y-daimeter.csv",row.names = FALSE)
write.csv(gene_anno_center,file="2020-1/siggene/Y-center.csv",row.names = FALSE)


