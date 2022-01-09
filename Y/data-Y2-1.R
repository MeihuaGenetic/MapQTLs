setwd("D:/Rdocument/mei-flower/Y")

LR_value_ML<-function(a,geno,pheno){
  L_colour <- c()
  pheno1<-as.matrix(as.numeric(pheno[,a]),ncol=1)
  for(i in 1:dim(geno)[2]){
    SNP <- (geno)[,i]
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

LR_size <-LR_value_ML(2,dat_Y2_st$geno,dat_Y2_st$pheno)
LR_pedicels <-LR_value_ML(4,dat_Y2_st$geno,dat_Y2_st$pheno)


LR_ML_permu<-function(a,geno,pheno,marker3){
  L1<-c()
  for(j in 1:100){
    L_colour<-c()
    pheno2<-as.matrix(as.numeric(pheno[,a]),ncol=1)
    num<-sample(nrow( pheno2))
    pheno1<-as.matrix(pheno2[num],ncol=1)
    for(i in marker3){
      SNP <- (geno)[,i]
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


plot(LR_size[geno_type_st$marker2])
plot(LR_size[geno_type_st$marker3])
LR_size_permu2 <-LR_ML_permu(2,dat_Y2_st$geno,dat_Y2_st$pheno,geno_type_st$marker2)# 6.827272
LR_size_permu3 <-LR_ML_permu(2,dat_Y2_st$geno,dat_Y2_st$pheno,geno_type_st$marker3)#6.072412


plot(LR_pedicels[geno_type_st$marker2])
plot(LR_pedicels[geno_type_st$marker3])
LR_pedicels_permu2 <-LR_ML_permu(4,dat_Y2_st$geno,dat_Y2_st$pheno,geno_type_st$marker2) #6.814851
LR_pedicels_permu3 <-LR_ML_permu(4,dat_Y2_st$geno,dat_Y2_st$pheno,geno_type_st$marker3)#6.041480


