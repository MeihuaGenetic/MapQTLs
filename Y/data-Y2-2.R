setwd("D:/Rdocument/mei-flower/Y")

dat_Y2.load <- function(geno_table=Y_marker_st$marker,pheno_table="Y-2.txt" ){
  dat_Y2<-list(TIME=NULL,
               geno=NULL,
               pheno=NULL,
               trait=1)
  A_F<-Y_marker_st$marker
  A_F<-t(A_F)
  A_F<-rbind(A_F[-1:-2,])
  
  B_F<-read.table("Y-2.txt",header=TRUE)
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
  
  dat_Y2$geno <- AA1_F.1
  dat_Y2$pheno <- BB1_F.1
  return(dat_Y2)
}
dat_Y2_st <- dat_Y2.load(geno_table=Y_marker_st$marker,pheno_table="Y-2.txt")

#############################################################################################


hist(dat_Y2_st$pheno[,1])
hist(dat_Y2_st$pheno[,2])
hist(dat_Y2_st$pheno[,3])
hist(dat_Y2_st$pheno[,4])


LR_mul_value<-function(a,geno,pheno,d){
  L_colour<-c()
  pheno11<-as.matrix(as.numeric(pheno[,a]),ncol=1)
  if(d==2){
    group<-c(0,(max(pheno11)-min(pheno11))/2,max(pheno11))
    pheno1<-cut(pheno11,breaks=group,labels = 1:(length(group)-1))
    pheno1<-as.matrix(as.numeric(pheno1),ncol=1)
  }
  if(d==3){
    group<-c(0,(max(pheno11)-min(pheno11))/3,2*(max(pheno11)-min(pheno11))/3,max(pheno11))
    pheno1<-cut(pheno11,breaks=group,labels = 1:(length(group)-1))
    pheno1<-as.matrix(as.numeric(pheno1),ncol=1)
  }
  if(d==4){
    group<-c(0,(max(pheno11)-min(pheno11))/4,(max(pheno11)-min(pheno11))/2,3*(max(pheno11)-min(pheno11))/4,max(pheno11))
    pheno1<-cut(pheno11,breaks=group,labels = 1:(length(group)-1))
    pheno1<-as.matrix(as.numeric(pheno1),ncol=1)
  }
  for(i in 1:dim(geno)[2]){
    SNP_F <- (geno)[,i]
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
  p
}

LR_pistils <-LR_mul_value(1,dat_Y2_st$geno,dat_Y2_st$pheno,2)
LR_petals <-LR_mul_value(3,dat_Y2_st$geno,dat_Y2_st$pheno,2)


LR_mul_permu<-function(a,geno,pheno,marker3,d){
  L1<-c()
  for(j in 1:100){
    L_colour<-c()
    pheno2 <- as.matrix(as.numeric(pheno[,a]),ncol=1)
    num<-sample(nrow( pheno2))
    pheno11<-as.matrix(pheno2[num],ncol=1)
    if(d==2){
      group<-c(0,(max(pheno11)-min(pheno11))/2,max(pheno11))
      pheno1<-cut(pheno11,breaks=group,labels = 1:(length(group)-1))
      pheno1<-as.matrix(as.numeric(pheno1),ncol=1)
    }
    if(d==3){
      group<-c(0,(max(pheno11)-min(pheno11))/3,2*(max(pheno11)-min(pheno11))/3,max(pheno11))
      pheno1<-cut(pheno11,breaks=group,labels = 1:(length(group)-1))
      pheno1<-as.matrix(as.numeric(pheno1),ncol=1)
    }
    if(d==4){
      group<-c(0,(max(pheno11)-min(pheno11))/4,(max(pheno11)-min(pheno11))/2,3*(max(pheno11)-min(pheno11))/4,max(pheno11))
      pheno1<-cut(pheno11,breaks=group,labels = 1:(length(group)-1))
      pheno1<-as.matrix(as.numeric(pheno1),ncol=1)
    }
    for(i in marker3){
      SNP_F <- (geno)[,i]
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

plot(LR_pistils[geno_type_st$marker2])
plot(LR_pistils[geno_type_st$marker3])
LR_pistils_permu2<-LR_mul_permu(1,dat_Y2_st$geno,dat_Y2_st$pheno,geno_type_st$marker2,2)#  7.131332 
LR_pistils_permu3<-LR_mul_permu(1,dat_Y2_st$geno,dat_Y2_st$pheno,geno_type_st$marker3,2)#5.101051 


plot(LR_petals[geno_type_st$marker2])
plot(LR_petals[geno_type_st$marker3])
LR_petals_permu2<-LR_mul_permu(3,dat_Y2_st$geno,dat_Y2_st$pheno,geno_type_st$marker2,2)# 7.029034 
LR_petals_permu3<-LR_mul_permu(3,dat_Y2_st$geno,dat_Y2_st$pheno,geno_type_st$marker3,2)#  4.920261


