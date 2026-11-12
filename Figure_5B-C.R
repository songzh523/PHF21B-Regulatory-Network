library(ggplot2)
library(ggpubr)
library(pheatmap)
library(psych)
library(dplyr)
library(corrplot)
library(openxlsx)

subloc=#significant genes from Jia L et al. Science. 2026 (science.adw1803_data_s7.xlsx) and our result (markergene_list.xlsx); columns must contain "cluster","gene","avg_log2FC"
##correlation between our RG subtypes and RG subtypes from Jia L et al. Science. 2026-svz and vz####
oursub<-read.xlsx(paste0(subloc,"markergene_list.xlsx"))
scisub<-read.xlsx(paste0(subloc,"science.adw1803_data_s7.xlsx"))

ourinp0<-oursub[,c("cluster","gene","avg_log2FC")]
sciinp0<-scisub[,c("cluster","gene","avg_log2FC")]

result<-matrix(NA,nc=3,nr=length(unique(ourinp0[,1]))*length(unique(ourinp0[,1])))
colnames(result)<-c("nam","pvalue","cor")
ids<-0
for(i in 1:length(unique(ourinp0[,1]))){ 
  for(j in 1:length(unique(sciinp0[,1]))){
    ids<-ids+1
    tmp1<-ourinp0[ourinp0[,1]==unique(ourinp0[,1])[i],]
    colnames(tmp1)[2]="gene"
    
    tmp2<-sciinp0[sciinp0[,1]==unique(sciinp0[,1])[j],]
    tmp<-merge(tmp1,tmp2,by="gene")
    
    cor1<-cor.test(tmp[,3],tmp[,5],method="spearman") 
    result[ids,1]=paste0(unique(ourinp0[,1])[i],"-",unique(sciinp0[,1])[j])
    result[ids,2]=cor1$p.value
    result[ids,3]=cor1$estimate[[1]]
    
  }
}
result<-as.data.frame(na.omit(result))
result$our<-sub("-.*","",result$nam)
result$sci<-sub(".*-","",result$nam)
mats0<-dcast(result[,c(3:ncol(result))],our~sci,value.var="cor")
mats<-mats0[,-1]
mats<-as.data.frame(lapply(mats,as.numeric))
rownames(mats)<-mats0[,1]
pmats0<-dcast(result[,c(2,4:ncol(result))],our~sci,value.var="pvalue")
pmats<-pmats0[,-1]
pmats<-as.data.frame(lapply(pmats,as.numeric))
rownames(pmats)<-pmats0[,1]
if (!is.null(pmats)){
  ssmt <- pmats< 0.01
  pmats[ssmt] <-'**'
  smt <- pmats >0.01& pmats <0.05
  pmats[smt] <- '*'
  pmats[!ssmt&!smt]<- ''
} else {
  pmats <- F
}

pheatmap(mats,scale = "none",cluster_row = T, cluster_col = T, border=NA,
         display_numbers = pmats,fontsize_number = 12, number_color = "black",
         cellwidth = 20, cellheight =20)

##correlation between our RG subtypes and RG subtypes from science paper-just vz####
oursub<-read.xlsx(paste0(subloc,"markergene_list.xlsx"))
scisub<-read.xlsx(paste0(subloc,"science.adw1803_data_s5.xlsx"))

ourinp0<-oursub[,c("cluster","gene","avg_log2FC")]
sciinp0<-scisub[,c("cluster","gene","avg_log2FC")]
result<-matrix(NA,nc=3,nr=length(unique(ourinp0[,1]))*length(unique(ourinp0[,1])))
colnames(result)<-c("nam","pvalue","cor")
ids<-0
for(i in 1:length(unique(ourinp0[,1]))){ 
  for(j in 1:length(unique(sciinp0[,1]))){
    ids<-ids+1
    tmp1<-ourinp0[ourinp0[,1]==unique(ourinp0[,1])[i],]
    colnames(tmp1)[2]="gene"
    
    tmp2<-sciinp0[sciinp0[,1]==unique(sciinp0[,1])[j],]
    tmp<-merge(tmp1,tmp2,by="gene")
    
    cor1<-cor.test(tmp[,3],tmp[,5],method="spearman") 
    result[ids,1]=paste0(unique(ourinp0[,1])[i],"-",unique(sciinp0[,1])[j])
    result[ids,2]=cor1$p.value
    result[ids,3]=cor1$estimate[[1]]
    
  }
}
result<-as.data.frame(na.omit(result))
result$our<-sub("-.*","",result$nam)
result$sci<-sub(".*-","",result$nam)
mats0<-dcast(result[,c(3:ncol(result))],our~sci,value.var="cor")
mats<-mats0[,-1]
mats<-as.data.frame(lapply(mats,as.numeric))
rownames(mats)<-mats0[,1]
pmats0<-dcast(result[,c(2,4:ncol(result))],our~sci,value.var="pvalue")
pmats<-pmats0[,-1]
pmats<-as.data.frame(lapply(pmats,as.numeric))
rownames(pmats)<-pmats0[,1]
if (!is.null(pmats)){
  ssmt <- pmats< 0.01
  pmats[ssmt] <-'**'
  smt <- pmats >0.01& pmats <0.05
  pmats[smt] <- '*'
  pmats[!ssmt&!smt]<- ''
} else {
  pmats <- F
}

pheatmap(mats,scale = "none",cluster_row = T, cluster_col = T, border=NA,
         display_numbers = pmats,fontsize_number = 12, number_color = "black",
         cellwidth = 20, cellheight =20)
