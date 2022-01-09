setwd("D:/Rdocument/mei-flower/Y")

dat_Y1.load <- function(geno_table=Y_marker_st$marker,pheno_table="Y-1.txt" ){
  dat_Y1<-list(TIME=NULL,
               geno=NULL,
               pheno=NULL,
               trait=1)
  A_F<-Y_marker_st$marker
  A_F<-t(A_F)
  A_F<-rbind(A_F[-1:-2,])
  
  B_F<-read.table("Y-1.txt",header=TRUE)
  B_F<-as.matrix(B_F)
  A1_F<-A_F[,1]
  B1_F<-B_F[,1]
  AB_F<- table(c(A1_F,B1_F))
  raw.id <- as.numeric(names(AB_F))
  C_F<-which(AB_F==2)
  ID_F<- sort(raw.id[C_F])
  
  nn_F <- length(ID_F)
  AA_F<- c()
  BB_F<- c()
  for(i in 1:nn_F){
    index2_F <- which(A_F[,1]==ID_F[i])
    AA_F <- rbind(AA_F,A_F[index2_F,])
    index3_F <- which(B_F[,1]==ID_F[i])
    BB_F <- rbind(BB_F,B_F[index3_F,])
  }   
  
  rownames(AA_F) <- AA_F[,1]
  AA1_F.1 <- AA_F[,-1]
  
  rownames(BB_F) <- BB_F[,1]
  BB1_F.1 <- BB_F[,-1]
  
  dat_Y1$geno <- AA1_F.1
  dat_Y1$pheno <- BB1_F.1
  return(dat_Y1)
}
dat_Y1_st <- dat_Y1.load(geno_table=Y_marker_st$marker,pheno_table="Y-1.txt")

geno_type.load<-function(geno=dat_Y1_st$geno){
  geno_type<-list(marker2=NULL,
                  marker3=NULL)
  
  marker_F_2<-c()
  marker_F_3<-c()
  for(i in 1:dim(dat_Y1_st$geno)[2]){
    SNP_F <- (dat_Y1_st$geno)[,i]
    SNP_F<-as.character(SNP_F)
    snp.type_F <- names(table(SNP_F))
    miss.type_F <- which(snp.type_F=="--")
    if(length(miss.type_F)>0){
      snp.type_F <- snp.type_F[-miss.type_F]
    }else{
      snp.type_F <- snp.type_F
    }
    if (length(snp.type_F)==2){
      marker_F_2<-c(marker_F_2,i)
    }else{
      marker_F_3<-c(marker_F_3,i)
    }
  }
  geno_type$marker2<-marker_F_2
  geno_type$marker3<-marker_F_3
  return(geno_type)
  
}
geno_type_st<-geno_type.load(geno=dat_Y1_st$geno)


#####################################################################
require(foreign)
require(nnet)
require(reshape2)
library(rms)


LR_value<-function(a,geno=dat_Y1_st$geno,pheno=dat_Y1_st$pheno){
  L_colour<-c()
  pheno1<-as.matrix(as.numeric(dat_Y1_st$pheno[,a]),ncol=1)
  for(i in 1:dim(dat_Y1_st$geno)[2]){
    SNP_F <- (dat_Y1_st$geno)[,i]
    SNP_F <- as.character(SNP_F)
    snp.type_F <- names(table(SNP_F))
    
    miss.type_F <- which(snp.type_F=="--")
    if(length(miss.type_F)>0){
      snp.type_F <- snp.type_F[-miss.type_F]
    }else{
      snp.type_F <- snp.type_F
    }
    miss_snp_F<-which(SNP_F=="--")
    if(length( miss_snp_F)>0){
      pheno3 <- pheno1[-miss_snp_F]
    }else{
      pheno3 <- pheno1
    }
    
    if(length( miss_snp_F)>0){
      SNP <- SNP_F[-miss_snp_F]
    }else{
      SNP <- SNP_F
    }
    
    SNP[which(SNP==snp.type_F[1])]<-c(1)
    SNP[which(SNP==snp.type_F[2])]<-c(2)
    SNP[which(SNP==snp.type_F[3])]<-c(3)
    SNP<-as.matrix(as.numeric(SNP),ncol=1)
    
    orm<-multinom(pheno3~SNP,trace = FALSE)
    orm2<-multinom(pheno3~1,trace = FALSE)
    
    LR<-orm2$deviance-orm$deviance
    
    snp.index <- length(snp.type_F)
    if (snp.index==3){
      df_1<-2
    }else{
      df_1<-1
    }
    L_1<-pchisq(LR,df=df_1,lower.tail=FALSE)
    L_colour<-c(L_colour,L_1)
    p<--1*log(L_colour)
  }
  p
}

LR_colour <-LR_value(1,dat_Y1_st$geno,dat_Y1_st$pheno)
LR_bud    <-LR_value(2,dat_Y1_st$geno,dat_Y1_st$pheno)
LR_center <-LR_value(3,dat_Y1_st$geno,dat_Y1_st$pheno)
LR_shape  <-LR_value(4,dat_Y1_st$geno,dat_Y1_st$pheno)
LR_period <-LR_value(5,dat_Y1_st$geno,dat_Y1_st$pheno)


LR_permu<-function(a,geno=dat_Y1_st$geno,pheno=dat_Y1_st$pheno,ID=geno_type_st$marker3){
  L1<-c()
  for(j in 1:100){
    L_colour<-c()
    pheno2<-as.matrix(as.numeric(dat_Y1_st$pheno[,a]),ncol=1)
    num<-sample(nrow( pheno2))
    pheno1<-as.matrix(pheno2[num],ncol=1)
    for(i in geno_type_st$marker3){
      SNP_F <- (dat_Y1_st$geno)[,i]
      SNP_F <- as.character(SNP_F)
      snp.type_F <- names(table(SNP_F))
      
      miss.type_F <- which(snp.type_F=="--")
      if(length(miss.type_F)>0){
        snp.type_F <- snp.type_F[-miss.type_F]
      }else{
        snp.type_F <- snp.type_F
      }
      miss_snp_F<-which(SNP_F=="--")
      if(length( miss_snp_F)>0){
        pheno3 <- pheno1[-miss_snp_F]
      }else{
        pheno3 <- pheno1
      }
      
      if(length( miss_snp_F)>0){
        SNP <- SNP_F[-miss_snp_F]
      }else{
        SNP <- SNP_F
      }
      
      SNP[which(SNP==snp.type_F[1])]<-c(1)
      SNP[which(SNP==snp.type_F[2])]<-c(2)
      SNP[which(SNP==snp.type_F[3])]<-c(3)
      SNP<-as.matrix(as.numeric(SNP),ncol=1)
      
      orm<-multinom(pheno3~SNP,trace = FALSE)
      orm2<-multinom(pheno3~1,trace = FALSE)
      
      LR<-orm2$deviance-orm$deviance
      snp.index <- length(snp.type_F)
      if (snp.index==3)
        df_1<-2
      else
        df_1<-1
      L_1<-pchisq(LR,df=df_1,lower.tail=FALSE)
      L_colour<-c(L_colour,L_1)
      p<--1*log(L_colour)
    }
    L_max<-max(p)
    L1<-c(L1,L_max)
    cat("num=",j,"LR=",L_max,"\n")
  }
  L1
}

plot(LR_colour[geno_type_st$marker2])
plot(LR_colour[geno_type_st$marker3])
LR_colour_permu2<-LR_permu(1,dat_Y1_st$geno,dat_Y1_st$pheno,geno_type_st$marker2)#5.897320
LR_colour_permu3<-LR_permu(1,dat_Y1_st$geno,dat_Y1_st$pheno,geno_type_st$marker3)#  6.181876


plot(LR_bud[geno_type_st$marker2])
plot(LR_bud[geno_type_st$marker3])
LR_bud_permu2<-LR_permu(2,dat_Y1_st$geno,dat_Y1_st$pheno,geno_type_st$marker2)#4.757960
LR_bud_permu3<-LR_permu(2,dat_Y1_st$geno,dat_Y1_st$pheno,geno_type_st$marker3)# 4.672774


plot(LR_center[geno_type_st$marker2])
plot(LR_center[geno_type_st$marker3])
LR_center_permu2<-LR_permu(3,dat_Y1_st$geno,dat_Y1_st$pheno,geno_type_st$marker2)# 5.771643
LR_center_permu3<-LR_permu(3,dat_Y1_st$geno,dat_Y1_st$pheno,geno_type_st$marker3)# 5.997276


plot(LR_center[geno_type_st$marker2])
plot(LR_center[geno_type_st$marker3])
LR_shape_permu2<-LR_permu(4,dat_Y1_st$geno,dat_Y1_st$pheno,geno_type_st$marker2)# 6.563671
LR_shape_permu3<-LR_permu(4,dat_Y1_st$geno,dat_Y1_st$pheno,geno_type_st$marker3)# 6.136940


plot(LR_period[geno_type_st$marker2])
plot(LR_period[geno_type_st$marker3])
LR_period_permu2<-LR_permu(5,dat_Y1_st$geno,dat_Y1_st$pheno,geno_type_st$marker2)# 7.227870
LR_period_permu3<-LR_permu(5,dat_Y1_st$geno,dat_Y1_st$pheno,geno_type_st$marker3) #7.177984 
