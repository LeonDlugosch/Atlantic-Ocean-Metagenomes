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
CheckDependencies = function(packages = NULL){
  inst.packages = as.data.frame(installed.packages())
  for(i in 1:length(packages)){
    if(all(grepl(packages[i], inst.packages$Package)) == TRUE){
      install.packages(packages[i])
    }
    library(packages[i], character.only = T)
  }
}
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
SilhouettePlot = function(clust = NULL, d = NULL, c.method = "ward.D2", d.method = "bray", max.k = 10, main = NA){
  CheckDependencies(packages = c("vegan", "cluster"))
  
  asw = as.vector(as.numeric(nrow(t(d))))
  sil.dist = as.matrix(vegdist(t(d), method = d.method))
  
  for(k in 2:max.k){
    sil = silhouette(cutree(clust, k = k), sil.dist)
    asw[k-1] = summary(sil)$avg.width
  }
  
  kbest = which.max(asw)+1 
  plot(x = 2:10,y = asw, type = "h", ylab = "Average silhouette width", xlab = "Number of clusters", main = main, las = 1)
  points(x = kbest, y = max(asw), pch = 16, col = "red")
  text(x = kbest, y = max(asw), paste("Optimal k =", kbest, sep = " "), pos = 4, col = "red")
}
GroupDistance = function(clust = NULL, d = NULL, k = NA, d.method = "bray"){
  CheckDependencies(packages = c("vegan", "cluster"))
  
  dist = as.matrix(vegdist(t(d), method = d.method))
  cut = cutree(clust, k = k)
  if(d.method == "bray"){ylab = "Bray-Curtis dissimilarity"}
  if(d.method == "euc"){ylab = "Euclidean distance"}
  for(i in 1:k){
    if(i == 1){
      sub.out.vec = NULL
      sub.in.vec = NULL
    }
    sub.in = dist[which(cut == k), which(cut == k)]
    for(j in 2:ncol(sub.in)){
      sub.in.vec = c(sub.in.vec, sub.in[j:ncol(sub.in),j-1])
    }
    
    sub.out = dist[which(cut != k), which(cut == k)]
    for(j in 1:ncol(sub.out)){
      sub.out.vec = c(sub.out.vec, sub.out[,j])
    }
  }
  if(d.method == "bray"){ylim = c(0,1)}
  if(d.method == "euc"){ylim = c(min(plot.dist$dist - 0.1*mean(plot.dist$dist)), max(plot.dist$dist - 0.1*mean(plot.dist$dist)))}
  plot.dist = rbind(data.frame(Group = "Ingroup", dist = sub.in.vec), data.frame(Group = "Outgroup", dist = sub.out.vec))
  krusk = kruskal.test(plot.dist$dist ~ plot.dist$Group)
  
  boxplot(plot.dist$dist ~ plot.dist$Group, ylab = ylab, axes = F, xlab = "",
          ylim = ylim, main = paste("Kruskal-Wallis chi²: ", round(krusk$statistic, 2), "\np-value: ", ifelse(krusk$p.value > 0.001, round(krusk$p.value, 3), "< 0.001")))
  axis(1, at = c(1,2), labels = c(paste("Within group\nn = ", table(plot.dist$Group)[1], sep = ""), paste("Outside group\nn = ", table(plot.dist$Group)[2], sep = "")))
  axis(2, las = 2)
  return(plot.dist)
}

ClosePlotDevice = function(){
  while(!is.null(dev.list())){dev.off()}
}
PlotLm = function(y = NULL, x = NULL, d = 1,
                  main = "", xlab = "Explanatory variable", ylab = "Response variable", type = "p", method = "lines", poly.col = "#66666688",
                  col = "black", bg = "black", l1.col = "red", l2.col = "black", pch = 21, cex = 1, las = 2, lty = 2, lwd = 2, zero.correction = F,
                  axes = TRUE, ylim = c(min(y)-(max(y)*0.1),max(y)+(max(y)*0.1)), xlim = NULL, mod.pos = "front", plot.result = NULL, add = FALSE, predict = 0.95){
  df = data.frame(x = x, y = y)
  df = df[complete.cases(df),]
  n = length(df$y)/length(unique(df$x))
  mod = lm(df$y ~ poly(df$x, d, raw = T))
  s.mod = summary(mod) 
  r = s.mod$adj.r.squared
  
  f = summary(mod)$fstatistic
  p = pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) = NULL
  
  if(p < 0.001){p = "< 0.001"}else{p = paste(" = ",round(p,3))}
  fit.mod = as.data.frame(predict(mod, interval = "conf", level = predict))
  fit.mod$f = x
  fit.mod = fit.mod[order(fit.mod$f),]
  if(add == FALSE){
    lm.plot = plot(y=y, x=x, type = "n", xlab = xlab, ylab = ylab, pch = pch, bg = bg, cex = cex, col = col, 
                   main = paste(main), axes = axes, ylim = ylim, xlim = xlim, las = las)
  }
  
  if(mod.pos == "front"){points(y=df$y, x=df$x, pch = pch, bg = bg, cex = cex, col = col)}
  if(is.numeric(predict)){
    if(method == "lines"){
      lines(fit.mod[,"lwr"] ~ fit.mod$f, col = l2.col, lty = 3, lwd = 2)
      lines(fit.mod[,"upr"] ~ fit.mod$f, col = l2.col, lty = 3, lwd = 2)
    }
    if(method == "poly"){
      polygon(x = c(fit.mod$f, rev(fit.mod$f)), c(fit.mod[,"upr"], rev(fit.mod[,"lwr"])), col = poly.col, border = NA)
    }
  }
  if(zero.correction == T){
    fit.mod[,"fit"][fit.mod[,"fit"] < 0] = 0
  }
  lines(fit.mod[,"fit"] ~ fit.mod$f, col = l1.col, lty = lty, lwd = 3)
  if(mod.pos == "back"){points(y=df$y, x=df$x, pch = pch, bg = bg, cex = cex, col = col)}
  if(mod.pos == "none"){}
  if(!is.null(plot.result)){
    text = paste(paste("adj.R2 = ", round(r,3),", ", sep = ""), paste("p.val ", p, sep = ""))
    legend(plot.result, bty = "n", legend = text)
    
  }
}
Normalize1 = function(df){
  if(class(df) == "data.frame" |class(df) == "matrix"){
    for (i in 1:ncol(df)){
      if (i == 1) {df2 = df}
      df2[,i] = (df[,i]-min(df[,i]))/(max(df[,i])-min(df[,i]))
    }
    return(df2)}
  if(class(df) == "numeric"){
    df2 = (df-min(df))/(max(df)-min(df))
    return(df2)
  }
}
Normalize100 = function(df){
  for (i in 1:ncol(df)){
    SUMS = colSums(df,na.rm = T)
    if (i == 1) {df2 = df}
    if(SUMS[i] == 0){
      df2[,i] = 0
    }else{df2[,i] = (df[,i]/SUMS[i])}
  }
  return(df2*100)  
}
cpm = function(data = NULL, gene_length = NULL){
  data = data/gene_length 
  sf = colSums(data)/1000000
  for (i in 1:ncol(data)){
    data[,i] = data[,i]/sf[i]
  }
  return(data)
}

#################################################################################################################
#                                                 Metadata                                                      #
#################################################################################################################
setwd("F:/MG_Analysis_v2/Data/")
Meta = read.csv2("ANT28_MG_Meta.csv")
setwd("F:/MG_Analysis_v2/Data/Supplement")
diversity = read.csv2("TableS3_Diversity_v3.csv", row.names = 1)
#################################################################################################################
#                                                 WOA data                                                      #
#################################################################################################################
#### Download nitrate World Ocean atlas data from: https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/nitrate/csv/all/1.00/woa18_all_n00mn01.csv.gz 
#### Download phsphate World Ocean atlas data from: https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/phosphate/csv/all/1.00/woa18_all_p00mn01.csv.gz

setwd("F:/MG_Analysis_v2/Data/WOA")
WOA.nox = read.csv("WOA_mean_nitrate_2018_annual_1deg.csv")
names(WOA.nox) = gsub("X", "", names(WOA.nox))
r.coordinates = data.frame(Longitude = round(Meta$Longitude)+0.5, Latitude = round(Meta$Latitude)+0.5, anu.mean.no3.surface = NA, anu.mean.no3.20m = NA)
for (i in 1:nrow(r.coordinates)){
  if (nrow(WOA.nox[which(WOA.nox$LATITUDE == r.coordinates$Latitude[i] & WOA.nox$LONGITUDE == r.coordinates$Longitude[i]),]) == 0){next}
  r.coordinates[i,3] = WOA.nox[which(WOA.nox$LATITUDE == r.coordinates$Latitude[i] & WOA.nox$LONGITUDE == r.coordinates$Longitude[i]),3]
  r.coordinates[i,4] = WOA.nox[which(WOA.nox$LATITUDE == r.coordinates$Latitude[i] & WOA.nox$LONGITUDE == r.coordinates$Longitude[i]),7]
}
Meta$Nitrate_annual_mean_20m = r.coordinates[,4]

WOA.po4 = read.csv("WOA_mean_phosphate_2018_annual_1deg.csv")
names(WOA.po4) = gsub("X", "", names(WOA.po4))
r.coordinates = data.frame(Longitude = round(Meta$Longitude)+0.5, Latitude = round(Meta$Latitude)+0.5, anu.mean.no3.surface = NA, anu.mean.no3.20m = NA)
for (i in 1:nrow(r.coordinates)){
  if (nrow(WOA.po4[which(WOA.po4$LATITUDE == r.coordinates$Latitude[i] & WOA.po4$LONGITUDE == r.coordinates$Longitude[i]),]) == 0){next}
  r.coordinates[i,3] = WOA.po4[which(WOA.po4$LATITUDE == r.coordinates$Latitude[i] & WOA.po4$LONGITUDE == r.coordinates$Longitude[i]),3]
  r.coordinates[i,4] = WOA.po4[which(WOA.po4$LATITUDE == r.coordinates$Latitude[i] & WOA.po4$LONGITUDE == r.coordinates$Longitude[i]),7]
}
Meta$Phosphate_annual_mean_20m = r.coordinates[,4]
Meta$NP_ratio = Meta$Nitrate_annual_mean_20m/Meta$Phosphate_annual_mean_20m
#################################################################################################################
#                           AOM data & generation of taxonomic and functional prfiles                           #
#################################################################################################################
setwd("")
setwd("F:/MG_Analysis_v2/Data/AOM_Mapped_Datasets")
AOM = read.csv2(file = "AOM_mapped_counts_v3.csv", header = T, stringsAsFactors = F)
AOM[,18:39] = cpm(data = AOM[,18:39], gene_length = AOM$Gene_Length)
data.red = as.data.frame(AOM[which(!is.na(AOM$Domain) & !is.na(AOM$Knr)),])
names(data.red)[18:39] = Meta$Station 
f.AOM = SummarizeDataset(data.red[,18:39], data.red$Knr)
t.AOM = SummarizeDataset(data.red[,18:39], paste(data.red$Domain, data.red$Phylum, data.red$Class, data.red$Order, data.red$Family, data.red$Genus, data.red$Species, sep = ";"))

#################################################################################################################
#                                       Cluster analysis/validation & NMDS                                      #
#################################################################################################################
tax.nmds = metaMDS(t(t.AOM), distance = "bray", k=2, trymax = 999, center=T)
tax.nmds$stress

KO.nmds = metaMDS(t(f.AOM), distance = "bray", k = 2, trymax = 999, center=T)
KO.nmds$stress

nr.nmds = metaMDS(t(data.red[,18:39]), distance = "bray", k=2, try = 999, center=T)
nr.nmds$stress

nr.clust = hclust(d = vegdist(t(data.red[,18:39]), method = "bray"), method = "ward.D2")
ko.clust = hclust(d = vegdist(t(f.AOM), method = "bray"), method = "ward.D2")
tax.clust = hclust(d = vegdist(t(t.AOM), method = "bray"), method = "ward.D2")

ClosePlotDevice()
Cairo(file = "Figure_S_Optimal_k.svg",
      type = "svg",
      width = 14,
      height = 14,
      bg = "white",
      units = "in",
      dpi = 96)
par(mfrow = c(3,3))

SilhouettePlot(clust = tax.clust, d = t.AOM, max.k = 10, c.method = "ward.D2", d.method = "bray", main = "Tax data")
GroupDistance(clust = tax.clust, d = t.AOM, k = 2, d.method = "bray")
plot(tax.clust, xlab = "", main = "Ward.D2 cluster dendrogram")

SilhouettePlot(clust = ko.clust, d = f.AOM, max.k = 10, c.method = "ward.D2", d.method = "bray", main = "KO data")
GroupDistance(clust = ko.clust, d = f.AOM, k = 3, d.method = "bray")
plot(ko.clust, xlab = "", main = "Ward.D2 cluster dendrogram")

SilhouettePlot(clust = nr.clust, d = data.red[,18:39], max.k = 10, c.method = "ward.D2", d.method = "bray", main = "nr data")
GroupDistance(clust = nr.clust, d = data.red[,18:39], k = 5, d.method = "bray")
plot(nr.clust, xlab = "", main = "Ward.D2 cluster dendrogram")
ClosePlotDevice()



Meta$NRcluster = cutree(nr.clust, k = 5)
Meta$KOcluster = cutree(ko.clust, k = 3)
Meta$TAXcluster = cutree(tax.clust, k = 2)
### renaming of KO clusters (ordered 1 = cold, 2 =  yellow, 3 = warm)
Meta$KOcluster = ifelse(Meta$KOcluster == 1, 2, ifelse(Meta$KOcluster == 2, 1, Meta$KOcluster))

write.csv2(data.frame(Station = Meta$Station, NRcluster = Meta$NRcluster, TAXcluster = Meta$TAXcluster, KOcluster = Meta$KOcluster),
           file = "AOM_cluster.csv", row.names = F)

#################################################################################################################
#                                  Data interpolation and random forest models                                  #
#################################################################################################################
library(randomForest)
RF_meta = Meta[,c(9,12,13,19,26,27,79,80,2,20,16)]
#write.csv(file = "Supplement_EnvData_us_3.csv", x = RF_meta, row.names = F, quote = F)

#### linear interpolation of missing Environmental data
RF_meta$Chla_ugL = approx(RF_meta$Chla_ugL, n = 22)$y
RF_meta$POC = approx(RF_meta$POC, n = 22)$y
RF_meta$Nitrate_annual_mean_20m = approx(RF_meta$Nitrate_annual_mean_20m, n = 22)$y
RF_meta$Phosphate_annual_mean_20m = approx(RF_meta$Phosphate_annual_mean_20m, n = 22)$y

#### Random forest dataset generation
RF_meta$AbsLat = abs(RF_meta$Latitude)
names(RF_meta) = c("Latitude", "Temperature", "Salinity", "Chlorophyll_a", "POC", "TPN", "Nitrate", "Phosphate", "Province", "BiomassProduction", "BacteriaCellcount", "Abs._Latitude")
RF_meta$Temperature = round(RF_meta$Temperature,3)
RF_meta$Nitrate = round(RF_meta$Nitrate,3)
write.csv(file = "Supplement_EnvData_us.csv", x = RF_meta, row.names = F, quote = F)
ClosePlotDevice()
set.seed(1234)
model.tax = randomForest(
  formula =  Meta$TAXcluster ~ .,
  data = RF_meta[,c(2:5, 7:9)]
)
varImpPlot(model.tax) 
which.min(model.tax$mse)
getTree(model.tax)

model.ko = randomForest(
  formula =  Meta$KOcluster ~ .,
  data = RF_meta[,c(2:5, 7:9)]
)
varImpPlot(model.ko) 
which.min(model.ko$mse)

model.nr = randomForest(
  formula =  Meta$KOcluster ~ .,
  data = RF_meta[,c(2:5, 7:9)]
)
varImpPlot(model.nr)
which.min(model.nr$mse)
importance(model.nr, )
model.nr$importance
randomF = data.frame(IncNodePurityTAX = model.tax$importance, IncNodePurityKO = model.ko$importance, IncNodePurityNR = model.nr$importance)
names(randomF) = c("TAX", "KO", "NR")

#################################################################################################################
#                                               Figure 2 plot                                                   #
#################################################################################################################
ClosePlotDevice()
Cairo(file = "Figure_2_Richness_Diversity_NMDS_woLegend.svg",
      type = "svg",
      width = 21,
      height = 12,
      bg = "white",
      units = "in",
      dpi = 96)
layout(matrix(c(4,4,5,5,1,1,1,10,10,
                6,6,7,7,2,2,2,11,11,
                8,8,9,9,3,3,3,12,12), nrow = 3, ncol = 9, byrow= T))
par(mar=c(4.1, 4.1, 4.1, 2.1))
ordisurf(tax.nmds, Meta$Temperature_T, col = "darkgrey", main = "Taxonomic profile")
ordispider(tax.nmds, Meta$TAXcluster, lty = 2, col = "black", label = T)
abline(h=0, lty=3)
abline(v=0, lty=3)
points(tax.nmds, display = "sites", pch = 21, cex = 2, bg = as.character(Meta$Prov_Col), lwd = .5)
text(x = -.4, y = .2, paste("Stress = ", round(tax.nmds$stress,3), sep = ""))

ordisurf(KO.nmds, Meta$Temperature_T, col = "darkgrey", main = "KO profile")
ordispider(KO.nmds, Meta$KOcluster, lty = 2, col = "black", label = T)
abline(h=0, lty=3)
abline(v=0, lty=3)
points(KO.nmds, display = "sites", pch = 21, cex = 2, bg = as.character(Meta$Prov_Col), lwd = .5)
text(x = -.4, y = .2, paste("Stress = ", round(KO.nmds$stress,3), sep = ""))

ordisurf(nr.nmds, Meta$Temperature_T, col = "darkgrey", main = "Gene profile")
ordispider(nr.nmds, Meta$NRcluster, lty = 2, col = "black", label = T)
abline(h=0, lty=3)
abline(v=0, lty=3)
points(nr.nmds, display = "sites", pch = 21, cex = 2, bg = as.character(Meta$Prov_Col), lwd = .5)
text(x = -2, y = 0.7, paste("Stress = ", round(nr.nmds$stress,3), sep = ""))
######### NEW PLOTS WITH POLYGONS! WHOOOOP!
#### TAXA
TAX = data.frame(Lat = Meta$Latitude, TAX.R = Meta$Species_Richness, TAX.D = Meta$Species_TSH)
TAX[,c(2,3)] = Normalize1(TAX[,c(2,3)])
for(i in 1:2){
  if(i == 1){mod.central = list()}
  n = length(TAX[,i+1])/length(unique(TAX[,1]))
  mod = lm(TAX[,i+1] ~ poly(TAX[,1], 5, raw = T))
  s.mod = summary(mod) 
  r = s.mod$adj.r.squared
  if(i == 1){plot(y = TAX[,i+1], x = TAX[,1], ylab = "Normalized richness & diversity", xlab = "", type = "n", main = "", las = 1)}
  fit.mod = as.data.frame(predict(mod, interval = "conf", level = .95))
  fit.mod$f = TAX[,1]
  fit.mod = fit.mod[order(fit.mod$f),]
  mod.central[[i]] = fit.mod$fit
  polygon(x = c(fit.mod$f, rev(fit.mod$f)), c(fit.mod[,"upr"], rev(fit.mod[,"lwr"])), col = ifelse(i == 1, "#fe702244", "#8a83ad44"), border = NA)
}
lines(y = mod.central[[1]], x = TAX$Lat, col = "#fe7022", lwd = 2, lty = 2)
lines(y = mod.central[[2]], x = TAX$Lat, col = "#8a83ad", lwd = 2, lty = 2)
points(y = TAX$TAX.R, x = TAX$Lat, pch = 21, bg = "#fe702266", col = "#66666666")
points(y = TAX$TAX.D, x = TAX$Lat, pch = 22, bg = "#8a83ad66", col = "#66666666")
legend("bottomleft", legend = c("Richness","Shannon diversity"), pch = 15, col = c("#8a83ad66", "#fe702266"),lty = 2, pt.cex = 2.5, bty = "n")

TAX = data.frame(TMP = Meta$Temperature_T, TAX.R = Meta$Species_Richness, TAX.D = Meta$Species_TSH)
TAX[,c(2,3)] = Normalize1(TAX[,c(2,3)])
for(i in 1:2){
  if(i == 1){mod.central = list()}
  n = length(TAX[,i+1])/length(unique(TAX[,1]))
  TAX = TAX[order(TAX$TMP),]
  mod = lm(TAX[,i+1] ~ poly(TAX[,1], 5, raw = T))
  s.mod = summary(mod) 
  r = s.mod$adj.r.squared
  if(i == 1){plot(y = TAX[,i+1], x = TAX[,1], ylab = "", xlab = "", type = "n", main = "", las = 1)}
  fit.mod = as.data.frame(predict(mod, interval = "conf", level = .95))
  fit.mod$f = TAX[,1]
  fit.mod = fit.mod[order(fit.mod$f),]
  mod.central[[i]] = fit.mod$fit
  polygon(x = c(fit.mod$f, rev(fit.mod$f)), c(fit.mod[,"upr"], rev(fit.mod[,"lwr"])), col = ifelse(i == 1, "#fe702244", "#8a83ad44"), border = NA)
}
lines(y = mod.central[[1]], x = TAX$TMP, col = "#fe7022", lwd = 2, lty = 2)
lines(y = mod.central[[2]], x = TAX$TMP, col = "#8a83ad", lwd = 2, lty = 2)
points(y = TAX$TAX.R, x = TAX$TMP, pch = 21, bg = "#fe702266", col = "#66666666")
points(y = TAX$TAX.D, x = TAX$TMP, pch = 22, bg = "#8a83ad66", col = "#66666666")
#### KEGGs
KEGG = data.frame(Lat = Meta$Latitude, KEGG.R = Meta$KOG_Richness, KEGG.D = Meta$KOG_TSH)
KEGG[,c(2,3)] = Normalize1(KEGG[,c(2,3)])
for(i in 1:2){
  if(i == 1){mod.central = list()}
  n = length(KEGG[,i+1])/length(unique(KEGG[,1]))
  mod = lm(KEGG[,i+1] ~ poly(KEGG[,1], 5, raw = T))
  s.mod = summary(mod) 
  r = s.mod$adj.r.squared
  if(i == 1){plot(y = KEGG[,i+1], x = KEGG[,1], ylab = "Normalized richness & diversity", xlab = "", type = "n", main = "", las = 1)}
  fit.mod = as.data.frame(predict(mod, interval = "conf", level = .95))
  fit.mod$f = KEGG[,1]
  fit.mod = fit.mod[order(fit.mod$f),]
  mod.central[[i]] = fit.mod$fit
  polygon(x = c(fit.mod$f, rev(fit.mod$f)), c(fit.mod[,"upr"], rev(fit.mod[,"lwr"])), col = ifelse(i == 1, "#fe702244", "#8a83ad44"), border = NA)
}
lines(y = mod.central[[1]], x = KEGG$Lat, col = "#fe7022", lwd = 2, lty = 2)
lines(y = mod.central[[2]], x = KEGG$Lat, col = "#8a83ad", lwd = 2, lty = 2)
points(y = KEGG$KEGG.R, x = KEGG$Lat, pch = 21, bg = "#fe702266", col = "#66666666")
points(y = KEGG$KEGG.D, x = KEGG$Lat, pch = 22, bg = "#8a83ad66", col = "#66666666")

KEGG = data.frame(TMP = Meta$Temperature_T, KEGG.R = Meta$KOG_Richness, KEGG.D = Meta$KOG_TSH)
KEGG[,c(2,3)] = Normalize1(KEGG[,c(2,3)])
for(i in 1:2){
  if(i == 1){mod.central = list()}
  n = length(KEGG[,i+1])/length(unique(KEGG[,1]))
  KEGG = KEGG[order(KEGG$TMP),]
  mod = lm(KEGG[,i+1] ~ poly(KEGG[,1], 5, raw = T))
  s.mod = summary(mod) 
  r = s.mod$adj.r.squared
  if(i == 1){plot(y = KEGG[,i+1], x = KEGG[,1], ylab = "", xlab = "", type = "n", main = "", las = 1)}
  fit.mod = as.data.frame(predict(mod, interval = "conf", level = .95))
  fit.mod$f = KEGG[,1]
  fit.mod = fit.mod[order(fit.mod$f),]
  mod.central[[i]] = fit.mod$fit
  polygon(x = c(fit.mod$f, rev(fit.mod$f)), c(fit.mod[,"upr"], rev(fit.mod[,"lwr"])), col = ifelse(i == 1, "#fe702244", "#8a83ad44"), border = NA)
}
lines(y = mod.central[[1]], x = KEGG$TMP, col = "#fe7022", lwd = 2, lty = 2)
lines(y = mod.central[[2]], x = KEGG$TMP, col = "#8a83ad", lwd = 2, lty = 2)
points(y = KEGG$KEGG.R, x = KEGG$TMP, pch = 21, bg = "#fe702266", col = "#66666666")
points(y = KEGG$KEGG.D, x = KEGG$TMP, pch = 22, bg = "#8a83ad66", col = "#66666666")
#### NRs
NR = data.frame(Lat = Meta$Latitude, NR.R = Meta$nrAOM_seq_Richness, NR.D = Meta$nrAOM_seq_TSH)
NR[,c(2,3)] = Normalize1(NR[,c(2,3)])
for(i in 1:2){
  if(i == 1){mod.central = list()}
  n = length(NR[,i+1])/length(unique(NR[,1]))
  mod = lm(NR[,i+1] ~ poly(NR[,1], 5, raw = T))
  s.mod = summary(mod) 
  r = s.mod$adj.r.squared
  if(i == 1){plot(y = NR[,i+1], x = NR[,1], ylab = "Normalized richness & diversity", xlab = "Latitude [°N]", type = "n", main = "", las = 1)}
  fit.mod = as.data.frame(predict(mod, interval = "conf", level = .95))
  fit.mod$f = NR[,1]
  fit.mod = fit.mod[order(fit.mod$f),]
  mod.central[[i]] = fit.mod$fit
  polygon(x = c(fit.mod$f, rev(fit.mod$f)), c(fit.mod[,"upr"], rev(fit.mod[,"lwr"])), col = ifelse(i == 1, "#fe702244", "#8a83ad44"), border = NA)
}
lines(y = mod.central[[1]], x = NR$Lat, col = "#fe7022", lwd = 2, lty = 2)
lines(y = mod.central[[2]], x = NR$Lat, col = "#8a83ad", lwd = 2, lty = 2)
points(y = NR$NR.R, x = NR$Lat, pch = 21, bg = "#fe702266", col = "#66666666")
points(y = NR$NR.D, x = NR$Lat, pch = 22, bg = "#8a83ad66", col = "#66666666")

NR = data.frame(TMP = Meta$Temperature_T, NR.R = Meta$nrAOM_seq_Richness, NR.D = Meta$nrAOM_seq_TSH)
NR[,c(2,3)] = Normalize1(NR[,c(2,3)])
for(i in 1:2){
  if(i == 1){mod.central = list()}
  n = length(NR[,i+1])/length(unique(NR[,1]))
  NR = NR[order(NR$TMP),]
  mod = lm(NR[,i+1] ~ poly(NR[,1], 5, raw = T))
  s.mod = summary(mod) 
  r = s.mod$adj.r.squared
  if(i == 1){plot(y = NR[,i+1], x = NR[,1], ylab = "", xlab = "Temperatur [°C]", type = "n", main = "", las = 1)}
  fit.mod = as.data.frame(predict(mod, interval = "conf", level = .95))
  fit.mod$f = NR[,1]
  fit.mod = fit.mod[order(fit.mod$f),]
  mod.central[[i]] = fit.mod$fit
  polygon(x = c(fit.mod$f, rev(fit.mod$f)), c(fit.mod[,"upr"], rev(fit.mod[,"lwr"])), col = ifelse(i == 1, "#fe702244", "#8a83ad44"), border = NA)
}
lines(y = mod.central[[1]], x = NR$TMP, col = "#fe7022", lwd = 2, lty = 2)
lines(y = mod.central[[2]], x = NR$TMP, col = "#8a83ad", lwd = 2, lty = 2)
points(y = NR$NR.R, x = NR$TMP, pch = 21, bg = "#fe702266", col = "#66666666")
points(y = NR$NR.D, x = NR$TMP, pch = 22, bg = "#8a83ad66", col = "#66666666")


library(fmsb)
randomF.n = Normalize100(randomF)
randomFT = rbind(rep(30,7), rep(0,7), as.data.frame(t(randomF.n)))

par(mar=c(1, 1, 1, 1))
radarchart(randomFT[c(1:2,3),], axistype = 1, axislabcol = "black",
           cglty = 1, cglwd = 1.5, cglcol = "darkgrey",
           pcol = "#ca9530ff", pfcol = "#ca953088",
           seg = 3, caxislabels=paste(seq(0,30,10), "%", sep = ""))

radarchart(randomFT[c(1:2,4),], axistype = 1, axislabcol = "black",
           cglty = 1, cglwd = 1.5, cglcol = "darkgrey",
           pcol = "#3a4d7dff", pfcol = "#3a4d7d88",
           seg = 3, caxislabels=paste(seq(0,30,10), "%", sep = ""))

radarchart(randomFT[c(1:2,5),], axistype = 1, axislabcol = "black",
           cglty = 1, cglwd = 1.5, cglcol = "darkgrey",
           pcol = "#65a252ff", pfcol = "#65a25288",
           seg = 3, caxislabels=paste(seq(0,30,10), "%", sep = ""))
ClosePlotDevice()
#################################################################################################################
#                                               Province Legend                                                 #
#################################################################################################################
par(mfrow = c(1,1), mar=c(4.1, 4.1, 4.1, 2.1))
ANT28.Legend = rev(c("APLR", "ANTA", "SANT", "FKLD", "SATL", "WTRA", "NATR", "NAST", "NADR"))
Legend.colors = rev(c("#85b6bb", "#226874","#0a303b", "#78369f", "#d04131", "#d27120", "#f1b434", "#558a33", "#327250"))
par(xpd=T)
plot(0,0, type = "n", axes = F, ylab = "", xlab = "")
legend("center", pch = 21, pt.bg = Legend.colors, legend = ANT28.Legend, pt.cex = 2, cex = 1, horiz = F)
par(xpd=F)

#################################################################################################################
#                          Temperature and geographic distance between stations                                 #
#################################################################################################################
Meta.Distance = Meta[,c(1,8,9,12)]
for (i in 1:nrow(Meta.Distance)){
  if (i == 1){R = 6371
  t.df = NULL
  dist.df = NULL} 
  lon.origin = Meta.Distance[i,2]*pi/180
  lat.origin = Meta.Distance[i,3]*pi/180
  temp.origin = Meta.Distance[i,4]
  for (j in 1:nrow(Meta.Distance)){
    if (j == 1){dist.vec = NULL
    t.vec = NULL}
    lon = Meta.Distance[j,2]*pi/180
    lat = Meta.Distance[j,3]*pi/180
    temp = Meta.Distance[j,4]
    
    t = abs(temp.origin-temp)
    d = acos(sin(lat.origin)*sin(lat) + cos(lat.origin)*cos(lat) * cos(lon.origin-lon)) * R  
    
    t.vec = c(t.vec, t)
    dist.vec = c(dist.vec, d)
    if (j == nrow(Meta.Distance)){dist.vec[is.nan(dist.vec)] = 0
    dist.df = cbind(dist.df, dist.vec)
    t.df = cbind(t.df, t.vec)}
  }
  if (i == nrow(Meta.Distance)) {
    row.names(dist.df) = Meta.Distance$Station
    dist.df = as.data.frame(dist.df)
    names(dist.df) = Meta.Distance$Station
    dist.df[dist.df < 1] = 0
    
    row.names(t.df) = Meta.Distance$Station
    t.df = as.data.frame(t.df)
    names(t.df) = Meta.Distance$Station
  }
}
rm(Meta.Distance, d, dist.vec, i, j, lat, lon, lat.origin, lon.origin, R, t, t.vec, temp, temp.origin)

#################################################################################################################
#                                       Distance/temperature decay analysis                                     #
#################################################################################################################
nr.dist = as.matrix(vegdist(t(data.red[,18:39]), method = "bray"))
tax.dist = as.matrix(vegdist(t(t.AOM), method = "bray"))
KO.dist = as.matrix(vegdist(t(f.AOM), method = "bray"))

##### Calculation of geographic distance and Temperature difference between stations ##### 
Meta.Distance = Meta[,c(1,8,9,12)]
for (i in 1:nrow(Meta.Distance)){
  if (i == 1){R = 6371
  t.df = NULL
  dist.df = NULL} 
  lon.origin = Meta.Distance[i,2]*pi/180
  lat.origin = Meta.Distance[i,3]*pi/180
  temp.origin = Meta.Distance[i,4]
  for (j in 1:nrow(Meta.Distance)){
    if (j == 1){dist.vec = NULL
    t.vec = NULL}
    lon = Meta.Distance[j,2]*pi/180
    lat = Meta.Distance[j,3]*pi/180
    temp = Meta.Distance[j,4]
    
    t = abs(temp.origin-temp)
    d = acos(sin(lat.origin)*sin(lat) + cos(lat.origin)*cos(lat) * cos(lon.origin-lon)) * R  
    
    t.vec = c(t.vec, t)
    dist.vec = c(dist.vec, d)
    if (j == nrow(Meta.Distance)){dist.vec[is.nan(dist.vec)] = 0
    dist.df = cbind(dist.df, dist.vec)
    t.df = cbind(t.df, t.vec)}
  }
  if (i == nrow(Meta.Distance)) {
    row.names(dist.df) = Meta.Distance$Station
    dist.df = as.data.frame(dist.df)
    names(dist.df) = Meta.Distance$Station
    dist.df[dist.df < 1] = 0
    
    row.names(t.df) = Meta.Distance$Station
    t.df = as.data.frame(t.df)
    names(t.df) = Meta.Distance$Station
  }
}
rm(Meta.Distance, d, dist.vec, i, j, lat, lon, lat.origin, lon.origin, R, t, t.vec, temp, temp.origin)

for(i in 1:ncol(t.df)){
  if(i == 1){
    res = as.data.frame(matrix(ncol = 5, nrow = 0))
    names(res) = c("Geo_dist", "Temp_dist", "tax_dist", "KO_dist", "nr_dist")
  }
  tmp = data.frame(Geo_dist = dist.df[,i], Temp_dist = t.df[,i], tax_dist = tax.dist[,i], KO_dist = KO.dist[,i], nr_dist = nr.dist[,i])
  res = rbind(res, tmp)
}
res = unique(res)

Cairo(file="Distance_Temperature_decay.svg", 
      type="svg",
      bg = "white",
      units="in", 
      width=7, 
      height=10, 
      dpi=96)
par(mar=c(5.1,6.2,4.1,2.1))
par(mfrow = c(3,1))
plot(0, 0, type = "n", axes = F, xlab = "", ylab = "")
legend("center", legend = c("Genes", "KOs", "Taxon"), lty = 1, col = c("#64a252", "#394c7d", "#c9952f"), lwd = 2, horiz = T, bty = "n")
PlotLm(y  = res$tax_dist, x = res$Temp_dist, l1.col = "#c9952f", d = 2, main = "",
       ylab = "Bray-Curtis dissimilarity", xlab = "Temperature [°C]", lty = 1, 
       axes = F, xlim = c(0,28), ylim = c(0, 1),  cex = .8, mod.pos = "front",
       predict = F, pch = 16, col = "#c9952f88")
PlotLm(y  = res$KO_dist, x = res$Temp_dist, l1.col = "#394c7d", d = 2, main = "",
       ylab = "", xlab = "", add = T, lty = 1,
       axes = F, xlim = c(0,28), ylim = c(0, 1),  cex = .8, mod.pos = "front",
       predict = F, pch = 16, col = "#394c7d88",)
PlotLm(y  = res$nr_dist, x = res$Temp_dist, l1.col = "#64a252", d = 2, main = "",
       ylab = "", xlab = "", add = T, lty = 1,
       axes = F, xlim = c(0,28), ylim = c(0, 1),  cex = .8, mod.pos = "front",
       predict = F, pch = 16, col = "#64a25288",)
axis(1, at = c(seq(from = 0, to = 25, by = 5)), labels = seq(from = 0, to = 25, by = 5))
axis(2, at = c(0, .2, .4, .6, .8, 1), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), las= 2)
box()

PlotLm(y  = res$tax_dist, x = res$Geo_dist, l1.col = "#c9952f", d = 2, main = "",
       ylab = "Bray-Curtis dissimilarity", xlab = "Distance [km x 103]", lty = 1, 
       axes = F, xlim = c(0,13000), ylim = c(0, 1),  cex = .8, mod.pos = "front",
       predict = F, pch = 16, col = "#c9952f88")
PlotLm(y  = res$KO_dist, x = res$Geo_dist, l1.col = "#394c7d", d = 2, main = "",
       ylab = "", xlab = "", add = T, lty = 1,
       axes = F, xlim = c(0,13000), ylim = c(0, 1),  cex = .8, mod.pos = "front",
       predict = F, pch = 16, col = "#394c7d88",)
PlotLm(y  = res$nr_dist, x = res$Geo_dist, l1.col = "#64a252", d = 2, main = "",
       ylab = "", xlab = "", add = T, lty = 1,
       axes = F, xlim = c(0,13000), ylim = c(0, 1),  cex = .8, mod.pos = "front",
       predict = F, pch = 16, col = "#64a25288",)
axis(1, at = c(seq(from = 0, to = 12500, by = 2500)), labels = c("0.0", "2.5", "5.0", "7.5", "10.0", "12.5"))
axis(2, at = c(0, .2, .4, .6, .8, 1), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), las= 2)
box()
ClosePlotDevice()
