#################################################################################################################
#                                     Cleaning Workspace & loading packages                                     #
#################################################################################################################
rm(list = ls(all=T))
library(DESeq2)
library(apeglm)
library(Cairo)
library(dplyr)
#################################################################################################################
#                                             Custom functions                                                  #
#################################################################################################################
SummarizeDataset = function(df, by, method = "sum"){
  print("Summarizing dataset...")
  S = names(df)
  for (i in 1:(length(S))){
    if (i == 1){
      print(paste(i, "of", length(S)))
      df.temp = data.frame(Class = by, Count = df[,i])
      names(df.temp) = c("Class", "Count")
      if(method == "sum"){df.s = summarise(group_by(df.temp, Class), Total = sum(Count))}
      if(method == "mean"){df.s = summarise(group_by(df.temp, Class), Mean = mean(Count))}
      names(df.s)[i+1] = S[i] 
    } 
    if (i > 1) {
      print(paste(i, "of", length(S)))
      df.temp = data.frame(Class = by, Count = df[,i])
      names(df.temp) = c("Class", "Count")
      if(method == "sum"){df.temp = summarise(group_by(df.temp, Class), Total = sum(Count))}
      if(method == "mean"){df.temp = summarise(group_by(df.temp, Class), Mean = mean(Count))}
      df.s = full_join(df.s, df.temp, by = "Class", .keep_all = T)
      names(df.s)[i+1] = S[i]
    }
    if (i == length(S)) {
      df.s[is.na(df.s)] = 0
      print(paste("Done! :-)"))
      df.s = as.data.frame(df.s)
      df.s[,1] = as.character(df.s[,1])
      df.s[,1][is.na(df.s[,1])] = "unkown"
      rownames(df.s) = df.s[,1]
      df.s = df.s[,-1]
      return(df.s)
    }
  }
}
ExtractField = function(x, n, sep = "_"){
  return(unlist(lapply(strsplit(x, sep), "[", n)))
}
ClosePlotDevice = function(){
  while(!is.null(dev.list())){dev.off()}
}
#################################################################################################################
#                                         Differential abundances of KOs                                        #
#################################################################################################################
setwd("F:/MG_Analysis_v2/Data/AOM_Mapped_Datasets")
AOMr = read.csv2(file = "AOM_count_data.csv", header = T, stringsAsFactors = F)
AOMr = AOMr[which(!is.na(AOMr$Domain) & (!is.na(AOMr$Knr) | !is.na(AOMr$CAZy_EC))),]
FuncID = paste(AOMr$Knr, AOMr$A, AOMr$B, AOMr$C, AOMr$Gene, sep = ";")
AOMr_KO = SummarizeDataset(AOMr[,c(18:39)], by = FuncID)
AOMr_CAZyme = AOMr[which(!is.na(AOMr$CAZy_EC)),]
AOMr_CAZyme = SummarizeDataset(AOMr_CAZyme[,c(18:39)], by = AOMr_CAZyme$CAZy_EC)

names(AOMr_KO) = ExtractField(names(AOMr_KO), 2, "_")

Cluster = read.csv2("AOM_cluster.csv") #Clusterdata generated in previous scrip "AOM_ClusterAnalysis_NMDS_DiversityPlot_RandomForest.R"

Cluster$KOcluster = as.factor(Cluster$KOcluster)
Cluster$NRcluster = as.factor(Cluster$NRcluster)
Cluster$TAXcluster = as.factor(Cluster$TAXcluster)

x = as.character(c("1", "1", "2"))
y = as.character(c("2", "3", "3"))
combinations = data.frame(x = x, y = y)

dseqmeta = data.frame(KO_cluster = Cluster$KOcluster)
rownames(dseqmeta) = names(AOMr_KO)

alpha = 0.05
logfold = 2
dds = DESeqDataSetFromMatrix(countData = AOMr_KO, 
                             colData = dseqmeta, 
                             design = ~ KO_cluster)
dds = DESeq(dds)
dds$KO_cluster
dds$sizeFactor
dds$replaceable

for (i in 1:nrow(combinations)){
  if (i == 1){
    result.list = list()
  }
  A = combinations[i,1]
  B = combinations[i,2]
  ### Again be mindful of the "Group" thingy
  res = results(dds, contrast = c("KO_cluster", A, B), independentFiltering = TRUE, alpha = alpha, pAdjustMethod = "BH", parallel = T) ### values are adjusted by Benjamini-Hochberg p-adjustment
  cols = densCols(res$lfcSE, -log10(res$padj))
  # ClosePlotDevice()
  # Cairo(file=paste("KO_", gsub("_", "", A), "_", gsub("_", "", B),".svg", sep = ""),
  #       type="svg",
  #       bg = "white",
  #       units="in",
  #       width=12,
  #       height=8,
  #       dpi=96)
  plot(res$log2FoldChange, -log10(res$padj), bg = cols, pch = 21, panel.first = grid(),
       main = paste(A, B, sep = "-"), xlab ="Effect size: log2 (fold-change)", ylab = "-log10 (adjusted p-value)",
       cex = 0.8)
  abline(v = 0)
  abline(v = c(-2,2), col = "brown")
  abline(h = -log10(alpha), col = "brown")
  # ClosePlotDevice()
  
  df.sig = data.frame(seq = rownames(res),l2fc = res$log2FoldChange, padj = res$padj)
  #df.sig = df.sig[which(abs(df.sig$l2fc) >= logfold & df.sig$padj <= alpha),]
  result.list[[paste(A,B,sep = "-")]] = df.sig
}
for(i in 1:length(result.list)){
  result.list[[i]] = data.frame(Knr = ExtractField(result.list[[i]]$seq, 1, ";"),
                                A = ExtractField(result.list[[i]]$seq, 2, ";"),
                                B = ExtractField(result.list[[i]]$seq, 3, ";"),
                                C = ExtractField(result.list[[i]]$seq, 4, ";"),
                                Gene = ExtractField(result.list[[i]]$seq, 5, ";"),
                                l2fc = result.list[[i]]$l2fc, 
                                padj = result.list[[i]]$padj)
}
setwd("C:/Users/icbmadmin/Desktop/AOM Code/Data")
pathways = read.csv2("AOM_differential_abundance_pathways.csv")
diff.ab_12 = result.list[[1]]
diff.ab_13 = result.list[[2]]
diff.ab_23 = result.list[[3]]

table(abs(diff.ab_12$l2fc) > 2 & diff.ab_12$padj <= 0.05)[2]/nrow(diff.ab_12)
table(abs(diff.ab_13$l2fc) > 2 & diff.ab_13$padj <= 0.05)[2]/nrow(diff.ab_12)
table(abs(diff.ab_23$l2fc) > 2 & diff.ab_23$padj <= 0.05)[2]/nrow(diff.ab_12)

p.levels = levels(as.factor(pathways$Pathways))
for (i in 1:length(p.levels)){
  if(i == 1){
    p.res = as.data.frame(matrix(ncol = 11, nrow = 0))
    names(p.res) = c("Class", "Pathway", "Function", "Knr", "Gene", "l2fc_12", "l2fc_12_p", "l2fc_13", "l2fc_13_p", "l2fc_23", "l2fc_23_p")
  }
  p = pathways[which(pathways$Pathways == p.levels[i]),]
  for(j in 1:nrow(p)){
    genes = strsplit(p$Gene.group[j], split = ",")
    for (k in 1:length(genes[[1]])){
      c1_c2 = result.list[[1]][which(grepl(genes[[1]][k], result.list[[1]]$Gene)),]
      c1_c3 = result.list[[2]][which(grepl(genes[[1]][k], result.list[[2]]$Gene)),]
      c2_c3 = result.list[[3]][which(grepl(genes[[1]][k], result.list[[3]]$Gene)),]
      if(nrow(c1_c2) > 0){
        tmp = data.frame(Class = NA, Pathway = NA, Knr = c1_c2$Knr, Gene = c1_c2$Gene)
        tmp$Class = p$Class[1]
        tmp$Pathway = p$Pathways[1]
        tmp$Function = p$Function[j]
        tmp$l2fc_12 = c1_c2$l2fc
        tmp$l2fc_12_p = c1_c2$padj
        tmp$l2fc_13 = c1_c3$l2fc
        tmp$l2fc_13_p = c1_c3$padj
        tmp$l2fc_23 = c2_c3$l2fc
        tmp$l2fc_23_p = c2_c3$padj
        p.res = rbind(p.res, tmp)
      }
    }
  }
  p.res = p.res[,c(1,2,5,3,4,6:11)]
}
p.res[,c(6:11)] = round(p.res[,c(6:11)], 4)
names(p.res)[3] = "Substrate"
write.csv2(file = "dseq2_differential_pathway_abundances_functional_clusters.csv", x = p.res, row.names = F)
p.res.s = p.res 

for (i in 1:length(p.levels)){
  if(i == 1){
    tmp.vec = list()
    bar.list = list()
  }
  tmp.vec[["1-2"]] = rep(0, 6)  
  tmp.vec[["1-3"]] = rep(0, 6)  
  tmp.vec[["2-3"]] = rep(0, 6)  
  tmp.data = p.res[which(p.res$Pathway == p.levels[i]),]  
  
  for(k in 1:3){
    c = ifelse(k == 1, 6, ifelse(k == 2, 8, 10))
    tmp.vec[[k]][1] = sum(tmp.data[,c] >= 4)
    tmp.vec[[k]][2] = sum(tmp.data[,c] < 4 & tmp.data[,c] >= 2)
    tmp.vec[[k]][3] = sum(tmp.data[,c] < 2 & tmp.data[,c] >= 0)
    tmp.vec[[k]][4] = sum(tmp.data[,c] > -2 & tmp.data[,c] <= 0)
    tmp.vec[[k]][5] = sum(tmp.data[,c] > -4 & tmp.data[,c] <= -2)
    tmp.vec[[k]][6] = sum(tmp.data[,c] <= -4)
    
    if(k == 3){
      bar.list[[p.levels[i]]] = data.frame(c1c2 = tmp.vec[[1]], c1c3 = tmp.vec[[2]], c2c3 = tmp.vec[[3]])
    }
  }
  
}

length(bar.list)
data.frame(id = names(bar.list), nr = 1:10)
bar.list = bar.list[c(3,7,1,4,5,2,6,8,9,10)]
save = bar.list
#################################################################################################################
#                                         Differential abundances of CAZymes                                    #
#################################################################################################################
dds = DESeqDataSetFromMatrix(countData = AOMr_CAZyme, 
                             colData = dseqmeta, 
                             design = ~ KO_cluster)
dds = DESeq(dds)
dds$KO_cluster
dds$sizeFactor
dds$replaceable

for (i in 1:nrow(combinations)){
  if (i == 1){
    result.list = list()
  }
  A = combinations[i,1]
  B = combinations[i,2]
  ### Again be mindful of the "Group" thingy
  res = results(dds, contrast = c("KO_cluster", A, B), independentFiltering = TRUE, alpha = alpha, pAdjustMethod = "BH", parallel = T) ### values are adjusted by Benjamini-Hochberg p-adjustment
  cols = densCols(res$lfcSE, -log10(res$padj))
  # ClosePlotDevice()
  # Cairo(file=paste("CAZyme_", gsub("_", "", A), "_", gsub("_", "", B),".svg", sep = ""),
  #       type="svg",
  #       bg = "white",
  #       units="in",
  #       width=12,
  #       height=8,
  #       dpi=96)
  plot(res$log2FoldChange, -log10(res$padj), bg = cols, pch = 21, panel.first = grid(),
       main = paste(A, B, sep = "-"), xlab ="Effect size: log2 (fold-change)", ylab = "-log10 (adjusted p-value)",
       cex = 0.8)
  abline(v = 0)
  abline(v = c(-2,2), col = "brown")
  abline(h = -log10(alpha), col = "brown")
  # ClosePlotDevice()
  
  df.sig = data.frame(seq = rownames(res),l2fc = res$log2FoldChange, padj = res$padj)
  #df.sig = df.sig[which(abs(df.sig$l2fc) >= logfold & df.sig$padj <= alpha),]
  result.list[[paste(A,B,sep = "-")]] = df.sig
}

diff.ab_12 = result.list[[1]]
diff.ab_13 = result.list[[2]]
diff.ab_23 = result.list[[3]]
table(abs(diff.ab_12$l2fc) > 2 & diff.ab_12$padj <= 0.05)[2]/nrow(diff.ab_12)
table(abs(diff.ab_13$l2fc) > 2 & diff.ab_13$padj <= 0.05)[2]/nrow(diff.ab_12)
table(abs(diff.ab_23$l2fc) > 2 & diff.ab_23$padj <= 0.05)[2]/nrow(diff.ab_12)

c.levels = c("AA", "CBM", "CE", "GH", "GT", "PL")
for(i in 1:length(c.levels)){
  
  if(i == 1){
    c.res = as.data.frame(matrix(ncol = 8, nrow = 0))
    names(c.res) = c("Class", "CAZy_EC", "l2fc_12", "padj_12", "l2fc_13", "padj_13", "l2fc_23", "padj_23")
  }
  
  c1_c2 = result.list[[1]][which(grepl(c.levels[i], result.list[[1]]$seq)),]
  c1_c2$Class = c.levels[i]
  names(c1_c2)[1:3] = c("CAZy_EC", "l2fc_12", "padj_12")
  c1_c2 = c1_c2[,c(4,1:3)]
  
  c1_c3 = result.list[[2]][which(grepl(c.levels[i], result.list[[2]]$seq)),]
  c1_c3$Class = c.levels[i]
  names(c1_c3)[1:3] = c("CAZy_EC", "l2fc_13", "padj_13")
  
  c2_c3 = result.list[[3]][which(grepl(c.levels[i], result.list[[3]]$seq)),]
  c2_c3$Class = c.levels[i]
  names(c2_c3)[1:3] = c("CAZy_EC", "l2fc_23", "padj_23")
  tmp = cbind(c1_c2, c1_c3[,2:3], c2_c3[,2:3])
  c.res = rbind(c.res, tmp)
  
}  

for (i in 1:length(c.levels)){
  if(i == 1){
    tmp.vec = list()
    bar.list = list()
  }
  tmp.vec[["1-2"]] = rep(0, 6)  
  tmp.vec[["1-3"]] = rep(0, 6)  
  tmp.vec[["2-3"]] = rep(0, 6)  
  tmp.data = c.res[which(c.res$Class == c.levels[i]),]  
  
  for(k in 1:3){
    c = ifelse(k == 1, 3, ifelse(k == 2, 5, 7))
    tmp.vec[[k]][1] = sum(tmp.data[,c] >= 4)
    tmp.vec[[k]][2] = sum(tmp.data[,c] < 4 & tmp.data[,c] >= 2)
    tmp.vec[[k]][3] = sum(tmp.data[,c] < 2 & tmp.data[,c] >= 0)
    tmp.vec[[k]][4] = sum(tmp.data[,c] > -2 & tmp.data[,c] <= 0)
    tmp.vec[[k]][5] = sum(tmp.data[,c] > -4 & tmp.data[,c] <= -2)
    tmp.vec[[k]][6] = sum(tmp.data[,c] <= -4)
    
    if(k == 3){
      bar.list[[c.levels[i]]] = data.frame(c1c2 = tmp.vec[[1]], c1c3 = tmp.vec[[2]], c2c3 = tmp.vec[[3]])
    }
  }
}

bar.list = c(save, bar.list)
length(bar.list)
c2.col = c("#f6e096", "#f4cd45", "#dcab00") #temperate
c1.col = c("#a4dbf7", "#45b8f3", "#0080c3") #cold
c3.col = c("#eda08c", "#f15f3b", "#cf2a00") #warm

c1c2_cols = c(rev(c1.col), c2.col)
c1c3_cols = c(rev(c1.col), c3.col)
c2c3_cols = c(rev(c2.col), c3.col)

ClosePlotDevice()
Cairo(file="Differential_Pathway_CAZyme_Abundance.svg",
      type="svg",
      bg = "white",
      units="in",
      width=8,
      height=37,
      dpi=96)
par(mfrow = c(16,3), c(5.1, 4.1, 1.1, 2.1), cex.lab = .7)
for(i in 1:length(bar.list)){
  barplot(bar.list[[i]]$c1c2/sum(bar.list[[i]]$c1c2), col = c1c2_cols, names.arg = "", xlab = ifelse(i == 16, "log2 fold-change", ""),
          ylim = c(0,1), axes = F, ylab = names(bar.list)[i], main = "")
  axis(2, labels = c(0, 50, 100), at = c(0,.5,1), las = 2)
  barplot(bar.list[[i]]$c1c3/sum(bar.list[[i]]$c1c3), col = c1c3_cols, names.arg = "", xlab = ifelse(i == 16, "log2 fold-change", ""),
          ylim = c(0,1), axes = F, ylab = "", main = "")
  barplot(bar.list[[i]]$c2c3/sum(bar.list[[i]]$c2c3), col = c2c3_cols, names.arg = "", xlab = ifelse(i == 16, "log2 fold-change", ""),
          ylim = c(0,1), axes = F, ylab = "", main = "")

}
ClosePlotDevice()


