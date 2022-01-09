setwd("D:/Rdocument/mei-flower/L")

dat_L3.load <- function(geno_table=L_marker_st$marker,pheno_table="L3.txt" ){
  dat_L3<-list(TIME=NULL,
               geno=NULL,
               pheno=NULL,
               trait=1)
  A_F<-L_marker_st$marker
  A_F<-t(A_F)
  A_F<-rbind(A_F[-1:-2,])
  
  B_F<-read.table("L3.txt",header=TRUE)
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
  
  dat_L3$geno <- AA1_F.1
  dat_L3$pheno <- BB1_F.1
  return(dat_L3)
}
dat_L3_st <- dat_L3.load(geno_table=L_marker_st$marker,pheno_table="L3.txt")

#############################################################################################

LR_value_ML<-function(a,geno=dat_L2_st$geno,pheno=dat_L2_st$pheno){
  L_colour <- c()
  pheno1<-as.matrix(as.numeric(dat_L2_st$pheno[,a]),ncol=1)
  for(i in 1:dim(dat_L2_st$geno)[2]){
    SNP <- (dat_L2_st$geno)[,i]
    SNP<-as.character(SNP)
    snp.type <- names(table(SNP))
    
    miss.type <- which(snp.type =="--")
    if(length(miss.type)>0){
      snp.type <- snp.type[-miss.type]
    }else{
      snp.type <- snp.type
    }
    
    
    snp.index <- length(snp.type)
    if (snp.index==3){
      AA<-pheno1[which(SNP==snp.type[1])]
      Aa<-pheno1[which(SNP==snp.type[2])]
      aa<-pheno1[which(SNP==snp.type[3])]
      S<-c(AA,Aa,aa)
      AA_mean<-mean(AA)
      Aa_mean<-mean(Aa)
      aa_mean<-mean(aa)
      S_mean<-mean(S)
      sd_AA<-sum((AA-AA_mean)^2)
      sd_Aa<-sum((Aa-Aa_mean)^2)
      sd_aa<-sum((aa-aa_mean)^2)
      S_sd<-sqrt((sd_AA+sd_Aa+sd_aa)/length(S))
      sum_sd<-sqrt(sum((S-S_mean)^2/length(S)))
      S_sum<-sum(-1*log(dnorm(S,S_mean,sum_sd)))
      AA_sum<-sum(-1*log(dnorm(AA,AA_mean,S_sd)))
      Aa_sum<-sum(-1*log(dnorm(Aa,Aa_mean,S_sd)))
      aa_sum<-sum(-1*log(dnorm(aa,aa_mean,S_sd)))
      LR<-(-2*((AA_sum+Aa_sum+aa_sum)-S_sum))
    }else{
      AA<-pheno1[which(SNP==snp.type[1])]
      Aa<-pheno1[which(SNP==snp.type[2])]
      S<-c(AA,Aa)
      AA_mean<-mean(AA)
      Aa_mean<-mean(Aa)
      S_mean<-mean(S)
      sd_AA<-sum((AA-AA_mean)^2)
      sd_Aa<-sum((Aa-Aa_mean)^2)
      S_sd<-sqrt((sd_AA+sd_Aa)/length(S))
      sum_sd<-sqrt(sum((S-S_mean)^2/length(S)))
      S_sum<-sum(-1*log(dnorm(S,S_mean,sum_sd)))
      AA_sum<-sum(-1*log(dnorm(AA,AA_mean,S_sd)))
      Aa_sum<-sum(-1*log(dnorm(Aa,Aa_mean,S_sd)))
      LR<-(-2*((AA_sum+Aa_sum)-S_sum))
      
    }
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

LR_length <-LR_value_ML(1,dat_L3_st$geno,dat_L3_st$pheno)
LR_daimeter <-LR_value_ML(2,dat_L3_st$geno,dat_L3_st$pheno)


LR_ML_permu<-function(a,geno=dat_L2_st$geno,pheno=dat_L2_st$pheno,ID=geno_type_st$marker3){
  L1<-c()
  for(j in 1:100){
    L_colour<-c()
    pheno2<-as.matrix(as.numeric(dat_L2_st$pheno[,a]),ncol=1)
    num<-sample(nrow( pheno2))
    pheno1<-as.matrix(pheno2[num],ncol=1)
    for(i in geno_type_st$marker3){
      SNP <- (dat_L2_st$geno)[,i]
      SNP<-as.character(SNP)
      snp.type <- names(table(SNP))
      
      miss.type <- which(snp.type =="--")
      if(length(miss.type)>0){
        snp.type <- snp.type[-miss.type]
      }else{
        snp.type <- snp.type
      }
      
      
      snp.index <- length(snp.type)
      if (snp.index==3){
        AA<-pheno1[which(SNP==snp.type[1])]
        Aa<-pheno1[which(SNP==snp.type[2])]
        aa<-pheno1[which(SNP==snp.type[3])]
        S<-c(AA,Aa,aa)
        AA_mean<-mean(AA)
        Aa_mean<-mean(Aa)
        aa_mean<-mean(aa)
        S_mean<-mean(S)
        sd_AA<-sum((AA-AA_mean)^2)
        sd_Aa<-sum((Aa-Aa_mean)^2)
        sd_aa<-sum((aa-aa_mean)^2)
        S_sd<-sqrt((sd_AA+sd_Aa+sd_aa)/length(S))
        sum_sd<-sqrt(sum((S-S_mean)^2/length(S)))
        S_sum<-sum(-1*log(dnorm(S,S_mean,sum_sd)))
        AA_sum<-sum(-1*log(dnorm(AA,AA_mean,S_sd)))
        Aa_sum<-sum(-1*log(dnorm(Aa,Aa_mean,S_sd)))
        aa_sum<-sum(-1*log(dnorm(aa,aa_mean,S_sd)))
        LR<-(-2*((AA_sum+Aa_sum+aa_sum)-S_sum))
      }else{
        AA<-pheno1[which(SNP==snp.type[1])]
        Aa<-pheno1[which(SNP==snp.type[2])]
        S<-c(AA,Aa)
        AA_mean<-mean(AA)
        Aa_mean<-mean(Aa)
        S_mean<-mean(S)
        sd_AA<-sum((AA-AA_mean)^2)
        sd_Aa<-sum((Aa-Aa_mean)^2)
        S_sd<-sqrt((sd_AA+sd_Aa)/length(S))
        sum_sd<-sqrt(sum((S-S_mean)^2/length(S)))
        S_sum<-sum(-1*log(dnorm(S,S_mean,sum_sd)))
        AA_sum<-sum(-1*log(dnorm(AA,AA_mean,S_sd)))
        Aa_sum<-sum(-1*log(dnorm(Aa,Aa_mean,S_sd)))
        LR<-(-2*((AA_sum+Aa_sum)-S_sum))
      }
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


plot(LR_length[geno_type_st$marker2])
plot(LR_length[geno_type_st$marker3])
LR_length_permu2 <-LR_ML_permu(1,dat_L3_st$geno,dat_L3_st$pheno,geno_type_st$marker2)#
LR_length_permu3 <-LR_ML_permu(1,dat_L3_st$geno,dat_L3_st$pheno,geno_type_st$marker3)#

plot(LR_daimeter[geno_type_st$marker2])
plot(LR_daimeter[geno_type_st$marker3])
LR_daimeter_permu2 <-LR_ML_permu(2,dat_L3_st$geno,dat_L3_st$pheno,geno_type_st$marker2) #
LR_daimeter_permu3 <-LR_ML_permu(2,dat_L3_st$geno,dat_L3_st$pheno,geno_type_st$marker3)# 



