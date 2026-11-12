library(RcisTarget) #R3.6.3
library(DT)
library("webshot")

##motifs of 42 PHF21B binding genes####
phf21bbinding<- #shared between PHF21B binding genes and DEGs in oRG_Inhibitory
###motif info####
motifRankings <- importRankings("hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
motifRankings
rank_cols<-colnames(motifRankings@rankings)

data(motifAnnotations_hgnc)
motifAnnotations_hgnc[1:3,1:6]

#remove pseudo genes and lncRNA-42
removed <- setdiff(phf21bbinding, rank_cols)
genelists<-list("PHF21B"=phf21bbinding[phf21bbinding %in% rank_cols])

#calcAUC
motifs_AUC <- calcAUC(genelists, motifRankings)

###cisTarget analysis(one-step)####
motifEnrichmentTable_wGenes <- cisTarget(genelists, 
                                         motifRankings,
                                         motifAnnot=motifAnnotations_hgnc,
                                         highlightTFs = NULL,
                                         nesThreshold = 3,
                                         geneErnMethod = "iCisTarget", 
                                         nCores = 4,
                                         verbose = TRUE)

motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)

resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]    
###PLOT MOTIF####
dtable<-datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
                  escape = FALSE, # To show the logo
                  filter="top", options=list(pageLength=5))
html<-"dtable.html"
saveWidget(dtable,html)


