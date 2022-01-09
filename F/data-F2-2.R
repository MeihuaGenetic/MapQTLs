
geno=dat_F2_st$geno
pheno=dat_F2_st$pheno

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


LR_pistils <-LR_mul_value(1,dat_F2_st$geno,dat_F2_st$pheno,3)
LR_petals <-LR_mul_value(3,dat_F2_st$geno,dat_F2_st$pheno,2)


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
LR_pistils_permu2<-LR_mul_permu(1,dat_F2_st$geno,dat_F2_st$pheno,geno_type_st$marker2,3)#14.983331
LR_pistils_permu3<-LR_mul_permu(1,dat_F2_st$geno,dat_F2_st$pheno,geno_type_st$marker3,3)# 10.665697

plot(LR_petals[geno_type_st$marker2])
plot(LR_petals[geno_type_st$marker3])
LR_petals_permu2<-LR_mul_permu(3,dat_F2_st$geno,dat_F2_st$pheno,geno_type_st$marker2,2)# 12.227272
LR_petals_permu3<-LR_mul_permu(3,dat_F2_st$geno,dat_F2_st$pheno,geno_type_st$marker3,2)#7.612247




