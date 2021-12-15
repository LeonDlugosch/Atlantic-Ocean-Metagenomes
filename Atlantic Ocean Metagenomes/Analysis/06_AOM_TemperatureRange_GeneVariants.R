#################################################################################################################
#                                     Cleaning Workspace & loading packages                                     #
#################################################################################################################
rm(list = ls(all=T))
library(vegan)
library(Cairo)
library(dplyr)
#################################################################################################################
#                                             Custom functions                                                  #
#################################################################################################################
cpm = function(data = NULL, gene_length = NULL){
  data = data/gene_length 
  sf = colSums(data)/1000000
  for (i in 1:ncol(data)){
    data[,i] = data[,i]/sf[i]
  }
  return(data)
}
ExtractField = function(x, n, sep = "_"){
  return(unlist(lapply(strsplit(x, sep), "[", n)))
}
ClosePlotDevice = function(){
  while(!is.null(dev.list())){dev.off()}
}
#################################################################################################################
#                                                 Metadata                                                      #
#################################################################################################################
setwd("F:/MG_Analysis_v2/Data/")
Meta = read.csv2("ANT28_MG_Meta.csv")
#################################################################################################################
#                           AOM data & generation of taxonomic and functional prfiles                           #
#################################################################################################################
setwd("")
setwd("F:/MG_Analysis_v2/Data/AOM_Mapped_Datasets")
AOM = read.csv2(file = "AOM_mapped_counts_v3.csv", header = T, stringsAsFactors = F)
AOM[,18:39] = cpm(data = AOM[,18:39], gene_length = AOM$Gene_Length)
data.red = as.data.frame(AOM[which(!is.na(AOM$Domain) & !is.na(AOM$Knr)),])


Tax.table = as.data.frame(table(data.red$Species))
Tax.table = Tax.table[rev(order(Tax.table$Freq)),]
Tax.table$nr = 1:nrow(Tax.table)
Taxa = as.character(Tax.table[which(Tax.table$Freq >= 2500), 1])

for(i in 1:length(Taxa)){
  if(i == 1){
    temp.range.taxon = as.data.frame(matrix(nrow = 1, ncol = 11))
    names(temp.range.taxon) = c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species/Genome", "Range", "MaxAbundanceTemp", "MinTemp", "MaxTemp")
  }
  print(Taxa[i])
  tmp.taxon = data.red[which(data.red$Species == Taxa[i]),]
  tax.range = data.frame(Abundance = colSums(tmp.taxon[,18:39])/max(colSums(tmp.taxon[,18:39])), Temperature = Meta$Temperature_T)
  temp.range.taxon[i,1] = tmp.taxon$Domain[1]
  temp.range.taxon[i,2] = tmp.taxon$Phylum[1]
  temp.range.taxon[i,3] = tmp.taxon$Class[1]
  temp.range.taxon[i,4] = tmp.taxon$Order[1]
  temp.range.taxon[i,5] = tmp.taxon$Family[1]
  temp.range.taxon[i,6] = tmp.taxon$Genus[1]
  temp.range.taxon[i,7] = Taxa[i] 
  temp.range.taxon[i,8] = max(tax.range[which(tax.range$Abundance >= 0.15),2]) - min(tax.range[which(tax.range$Abundance >= 0.15),2])
  temp.range.taxon[i,9] = tax.range[which(max(tax.range$Abundance) == tax.range$Abundance), 2]  
  temp.range.taxon[i,10] = min(tax.range[which(tax.range$Abundance >= 0.15),2])
  temp.range.taxon[i,11] = max(tax.range[which(tax.range$Abundance >= 0.15),2])
}
temp.range.taxon$Domain[temp.range.taxon$Domain != "Eukaryota"] = "Prokaryotes"
kruskal.test(temp.range.taxon$Range ~ temp.range.taxon$Domain)

write.csv2(file = "Taxon_temperature_range.csv", temp.range.taxon, row.names = F)
taxa_ridge = c("Candidatus Pelagibacter ubique", "Candidatus Pelagibacter sp. HTCC7211", "Candidatus Pelagibacter sp. TMED128",
               "alpha proteobacterium HIMB5", "SAR116 cluster alpha proteobacterium HIMB100", "Rhodobacteraceae bacterium TMED111", "Planktomarina temperata",
               "SAR86 cluster bacterium SAR86A","SAR86 cluster bacterium SAR86B", "SAR86 cluster bacterium SAR86E","SAR92 bacterium BACL26 MAG-121220-bin70",
               "Gammaproteobacteria bacterium TMED112","Gammaproteobacteria bacterium TMED186","Gammaproteobacteria bacterium TMED236",
               "Cellvibrionales bacterium TMED122",
               "Synechococcus sp.", "Synechococcus sp. CC9605",
               "Prochlorococcus marinus", "Prochlorococcus sp.",
               "Euryarchaeota archaeon TMED99", "Euryarchaeota archaeon TMED103",
               "Bathycoccus prasinos", "Aureococcus anophagefferens","Micromonas commoda", "Micromonas pusilla", "Ostreococcus 'lucimarinus'")

setwd("F:/MG_Analysis_v2/Diversity/NatCom")
for (i in 1:length(taxa_ridge)){
  if(i == 1){
    Cairo(file = "Taxon_Range.svg",
          type = "svg",
          width = 12,
          height = 22,
          bg = "white",
          units = "in",
          dpi = 96)
    ridge_res = as.data.frame(matrix(ncol = 23, nrow = 1))
    names(ridge_res) = c("Taxon", Meta$Station)
    par(mfrow = c(9,3), mar = c(5.1,4,3.1,2.1), cex.main = .9)
  }
  tmp.taxon = data.red[which(data.red$Species == taxa_ridge[i]),]
  tax.range = data.frame(Abundance = colSums(tmp.taxon[,18:39])/max(colSums(tmp.taxon[,18:39])), Temperature = Meta$Temperature_T)
  ridge_res[i,c(2:23)] = colSums(tmp.taxon[,18:39])/max(colSums(tmp.taxon[,18:39]))
  ridge_res[i,1] = taxa_ridge[i]
  minT = min(tax.range[which(tax.range$Abundance >= 0.15),2])
  maxT = max(tax.range[which(tax.range$Abundance >= 0.15),2])
  # Add line on top
  tax.range = tax.range[order(tax.range$Temperature),]
  plot(x = tax.range$Temperature, tax.range$Abundance,
       col = "darkgrey", type = "l", lwd=1, ylim = c(-0.05,1),
       main = taxa_ridge[i],
       xlab = ifelse(i >= 24, "Temperature [°C]", ""),
       ylab = ifelse(i == 1 | i == 4 | i == 7 | i == 10 | i == 13 | i == 16 | i == 19 | i == 22 | i == 25, "Normalised abundance", ""),
       las = 1, axes = F)
  axis(1)
  axis(2, at = c(0, .5, 1), labels = c("0.0", "0.5", "1.0"))
  box()
  polygon( 
    c(min(tax.range$Temperature), tax.range$Temperature, max(tax.range$Temperature)), 
    c(0, tax.range$Abundance, 0), 
    col = "lightgrey", border=F)
  rect(minT, -0.05, maxT, 0, col = "#FF0000", border = NA)
  abline(h = .15, lty = 2)
  abline(v = tax.range[which(tax.range$Abundance == max(tax.range$Abundance)),2], lty = 1)
  if(i == length(taxa_ridge)){
    plot(0,0, axes = F, main = "", ylab = "", xlab = "", type = "n")
    legend("center", legend = c("Max. abundance", "15% of max. abundance", "Taxon temperature range"), bty = "n",
           xpd = T, pch = c(NA, NA, 15), lty = c(1,2,NA), col = c("black", "black", "red"))
    
    ClosePlotDevice()
  }
}

for(i in 1:length(Taxa)){
  if(i == 1){
    res = as.data.frame(matrix(nrow = 1, ncol = 6))
    names(res) = c("Taxon", "KO/Gene", "Number of variants", "Clusters", "mean temperature range", "sd")
    #temp.range.taxon = as.data.frame(matrix(nrow = 1, ncol = 2))
    c = 1
  }
  print(Taxa[i])
  tmp.taxon = data.red[which(data.red$Species == Taxa[i]),]
  tax.range = data.frame(Abundance = colSums(tmp.taxon[,18:39])/max(colSums(tmp.taxon[,18:39])), Temperature = Meta$Temperature_T)
  max(tax.range[which(tax.range$Abundance >= 0.15),2]) - min(tax.range[which(tax.range$Abundance >= 0.15),2])
  KO.table = as.data.frame(table(paste(tmp.taxon$Knr,tmp.taxon$Gene, sep = ";")))
  KO.table = KO.table[rev(order(KO.table$Freq)),]
  KO.table = KO.table[which(!grepl("rpo", KO.table$Var1)),]
  
  KO.table = as.character(KO.table[which(KO.table$Freq >= 10),1])
  KO.table = ExtractField(KO.table, 1, ";")
  
  for(j in 1:length(KO.table)){
    
    tmp.taxon.ko = tmp.taxon[which(tmp.taxon$Knr == KO.table[j]),]
    tmp.taxon.ko = tmp.taxon.ko[which(rowSums(tmp.taxon.ko[18:39]) > 0),]
    
    res[c,1] = Taxa[i]
    res[c,2] = paste(tmp.taxon.ko$Knr[1], tmp.taxon.ko$Gene[1], sep = ";")
    res[c,3] = nrow(tmp.taxon.ko)
    print(res[c,2])
    
    for(z in 1:nrow(tmp.taxon.ko)){
      if(z == 1){t.range = NULL}
      nr.range = data.frame(Abundance = as.numeric(tmp.taxon.ko[z,18:39]), Temperature = Meta$Temperature_T)
      t.range = c(t.range, max(nr.range[which(nr.range$Abundance >= 0.15),2]) - min(nr.range[which(nr.range$Abundance >= 0.15),2]))
    }
    
    res[c,5] = mean(t.range[t.range > 0])
    res[c,6] = sd(t.range[t.range > 0])
    c = c+1
  }
}
levels(as.factor(res$Taxon))
res.box = res[which(!grepl("unid.", res$Taxon) & !grepl("NA", res$Taxon, ignore.case = F)),]
levels(as.factor(res.box$Taxon))
boxplot(res.box$`mean temperature range` ~ res.box$Taxon, horizontal = T, las = 1, ylab = "", xlab = "Mean variant temperature range [°C]", col = "orange", pch = 16, cex =.2)

write.csv2(file = "Gene_variants_single_NR.csv", res, row.names = F)

for(i in 1:length(taxa_ridge)){
  if(i == 1){plot_range = res.box[which(res.box$Taxon == taxa_ridge[i]),]}else{plot_range = rbind(plot_range, res.box[which(res.box$Taxon == taxa_ridge[i]),])}
}
setwd("F:/Work/ANT28 Diversity/NatComms/Plots")
ClosePlotDevice()
Cairo(file = "Figure_5_TemperatureRange.svg",
      type = "svg",
      width = 18,
      height = 12,
      bg = "white",
      units = "in",
      dpi = 96)
par(mar = c(5,20,2,1), cex.axis = 1)
layout(matrix(c(1,1,1,2,2,
                1,1,1,3,3,
                1,1,1,3,3,
                1,1,1,3,3,
                1,1,1,3,3), nrow = 5, ncol = 5, byrow = T))
plot_range$Taxon = factor(plot_range$Taxon, levels = rev(taxa_ridge))
plot(0,0, xlim = c(0,28), ylim = (c(1,26)), xlab = "", ylab = "", axes = F, type = "n")
abline(h = c(5.5, 7.5, 11.5), col = "lightgrey")
abline(v = 8.59, col = "black")

boxplot(plot_range$`mean temperature range` ~ plot_range$Taxon,
        horizontal = T, las = 1, ylab = "Temperature range [°C]", xlab = "Mean variant temperature range [°C]", col = rev(c(rep("#fcc87c", 15), rep("#ec8356", 4), rep("#579797", 2), rep("#345969", 5))),
        pch = 16, cex =.2, ylim = c(0,28), add = T)
t_range = NULL
for(i in 1:length(taxa_ridge)){
  t_range = c(t_range, temp.range.taxon$Range[which(temp.range.taxon[,7] == taxa_ridge[i])])
}
points(y = 1:26, x = rev(t_range), cex = 2, pch = 23, bg = rev(c(rep("#fcc87c", 15),rep("#ec8356", 4), rep("#579797", 2),  rep("#345969", 5))))
#par(mfrow = c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1), cex.axis = 1)
plot(0,0, type = "n", axes = F, main = "", ylab = "", xlab = "")
legend("center", pch = c(23, 22, 22, 22, 22), pt.cex = c(2, 2.5, 2.5, 2.5, 2.5), cex = 2, pt.bg = c("white", "#fcc87c", "#ec8356", "#579797", "#345969"),
       legend = c("Taxon range", "Heterotrophic bacteria", "Cyanobacteria", "Archaea", "Eukaryota"), xpd = T)

Kref = read.delim("F:/MG_Analysis_v2/Data/KEGG Reference/K_Ref.txt", sep = "\t", header = T)
KEGG_Ranges = res[,c(2,5)]
KEGG_Ranges$`KO/Gene` = ExtractField(KEGG_Ranges$`KO/Gene`, 1, ";")
names(KEGG_Ranges)[1] = "Knr"

KEGG_Ranges = left_join(KEGG_Ranges, Kref)
#boxplot(KEGG_Ranges$`mean temperature range` ~ KEGG_Ranges$B)
KEGG_Ranges$B = trimws(KEGG_Ranges$B,which = "both")
levels(as.factor(KEGG_Ranges$B))
KEGG_Ranges$B[KEGG_Ranges$B == "Enzyme families"] = "Peptidases andy protein kinases"
KEGG_pathways = c("Translation", "Energy metabolism", "Carbohydrate metabolism", "Nucleotide metabolism","Amino acid metabolism",
                  "General Metabolism","Metabolism of cofactors and vitamins",
                  "Replication and repair", "Membrane transport", 
                  "Peptidases andy protein kinases")


for(i in 1:length(KEGG_pathways)){
  if(i == 1){Pathway_range = KEGG_Ranges[which(KEGG_Ranges$B == KEGG_pathways[i]),]}else{Pathway_range = rbind(Pathway_range,KEGG_Ranges[which(KEGG_Ranges$B == KEGG_pathways[i]),])}
}
Pathway_range$B = factor(Pathway_range$B, levels = rev(KEGG_pathways))
table(Pathway_range$B)
par(mar = c(5,1,2,20), cex.axis = 1)
plot(0,0, xlim = c(0,20), ylim = c(0.5,10.5), xlab = "Temperature range [°C]", ylab = "", axes = F, type = "n")
boxplot(Pathway_range$`mean temperature range` ~ Pathway_range$B, horizontal = T, las = 2, ylab = "", cex = NA, ylim = c(0,20),
        names = NA, add = T, axes = F)
abline(v = 8.59)
axis(1)
axis(4, labels = paste(rev(KEGG_pathways), " (n = ", table(Pathway_range$B), ")", sep = ""), at = 1:10, las = 2)
box()
ClosePlotDevice()

write.csv2(file = "Taxon_Range.csv", x = temp.range.taxon, row.names = F)
res$KO = ExtractField(res$`KO/Gene`, 1, ";")
res$Gene = ExtractField(res$`KO/Gene`, 2, ";")
names(res)[7] = "Knr"
res = left_join(res, Kref)
res = res[,c(1,7:11,3,5:6)]
write.csv2(file = "TaxKO_Range.csv", x = res, row.names = F)
