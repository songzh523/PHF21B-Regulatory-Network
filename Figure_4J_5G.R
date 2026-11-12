##disease gsva/gsea analysis-each region each celltype one by one-celltype####
library(ggplot2)
library(ggrepel)
library(liger)
library(stringr)
library(dplyr)
library(clusterProfiler)

resloc= #folders with gene fold changes between mutation and wild type organoids; all files with suffix "_sigmarkers_PandN_all-bimod.txt"
gmtloc= #folders with gene sets download from molecular signatures database (MSigDB): format txt files  

#1 make gmt files
gmtfiles<-list.files(path=gmtloc,pattern = "*.txt")
gmtfilesnam0<-gsub("\\..*","",gmtfiles)
geneselects<-c("ASTON_MAJOR_DEPRESSIVE_DISORDER_UP","KEGG_ALZHEIMERS_DISEASE","KEGG_HUNTINGTONS_DISEASE","ASTON_MAJOR_DEPRESSIVE_DISORDER_UP","HP_INTELLECTUAL_DISABILITY_SEVERE","HP_EPILEPTIC_ENCEPHALOPATHY","HP_ATTENTION_DEFICIT_HYPERACTIVITY_DISORDER","HP_AUTISM") #,33
gmtfilesnam<-gmtfilesnam0[gmtfilesnam0 %in% geneselects]
kepgmtfiles<-gmtfiles[gmtfilesnam0 %in% geneselects]
gs.lists<-list()

for(f in 1:length(kepgmtfiles)){
  gene_sets <- read.table(paste0(gmtloc,kepgmtfiles[f]))
  gene_names <- unlist(gene_sets)
  gs.lists[[f]]=gene_names
  names(gs.lists)[f]=gmtfilesnam[f]
}

listtomatrix<-function(uneven_list){
  for(i in 1:length(uneven_list)){
    if(i ==1){dat<-as.data.frame(uneven_list[[i]])
    dat$type<-names(uneven_list)[i]
    colnames(dat)[1]="gene"}else{
      tmp<-as.data.frame(uneven_list[[i]])
      tmp$type<-names(uneven_list)[i]
      colnames(tmp)[1]="gene"
      dat<-rbind(dat,tmp)
    }
  }
  return(dat)
}

genesets_mat0<-listtomatrix(gs.lists)
genesets_mat<-genesets_mat0[,c(2,1)]
colnames(genesets_mat)<-c("term","gene")

#2 read deg genes,and gsea analysis
files<-list.files(path = resloc,pattern = "*_sigmarkers_PandN_all-bimod.txt")
celltypes<-sub("-.*","",files)
for(c in 1:length(files)){
  markers<-read.table(paste0(resloc,files[c]))
  markers_mH<-markers
  markers_mH = markers_mH[order(markers_mH$avg_log2FC, decreasing = T),]
  genelist <- structure(markers_mH$avg_log2FC, names = rownames(markers_mH)) 
  res <- GSEA(genelist, TERM2GENE = genesets_mat, pvalueCutoff = 1,eps=0) 
  write.table(data.frame(res),paste0(resloc,celltypes[c],"_GSEA_PandN_all-bimod-disease"),sep="\t",row.names = F,quote = F)
  
}

#3 gsea analysis plot
fileuse<-list.files(path = resloc,pattern = "*_GSEA_PandN_all-bimod-disease$") 
for(f in 1:length(fileuse)){
  if(f==1){
    dats<-read.table(paste0(resloc1,fileuse[f]),sep="\t",header = T)
    dats$typ<-sub("_GSEA_PandN_all-bimod-disease","",fileuse[f])
  }else{
    tmp<-read.table(paste0(resloc1,fileuse[f]),sep="\t",header = T)
    if(nrow(tmp)>0){
      tmp$typ<-sub("_GSEA_PandN_all-bimod-disease","",fileuse[f])
      dats<-rbind(dats,tmp)
    }
  }
}
cmt0<-reshape2::dcast(dats[,c("Description","typ","NES")],Description~typ)
cmt<-cmt0[,-1]
rownames(cmt)<-cmt0[,1]
pmt0<-reshape2::dcast(dats[,c("Description","typ","pvalue")],Description~typ)
pmt<-pmt0[,-1]
rownames(pmt)<-pmt0[,1]
rownames(cmt)<-firstCap(gsub("_"," ",sub("GOBP_","",rownames(cmt))))
rownames(pmt)<-firstCap(gsub("_"," ",sub("GOBP_","",rownames(pmt))))
pmt1<-(-log10(pmt))

pheatmap(pmt1, display_numbers = matrix(ifelse(pmt <0.01, "**", ifelse(pmt <0.05,"*","")),nrow(pmt)),cluster_row=F,treeheight_col =0,color=colorRampPalette(c("#F3D5D3","#A32530"))(50),border="white")


