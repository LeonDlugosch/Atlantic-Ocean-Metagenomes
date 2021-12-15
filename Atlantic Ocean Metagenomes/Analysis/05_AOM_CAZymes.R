#################################################################################################################
#                                     Cleaning Workspace & loading packages                                     #
#################################################################################################################
rm(list = ls(all=T))
library(Cairo)
library(dplyr)
library(vegan)
library(RColorBrewer)
library(ape)
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
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
SumBool = function(x = NULL){
  bools = c(ifelse(length((x[which(x[,1] == F),2])) > 0,x[which(x[,1] == F),2],0), 
            ifelse(length((x[which(x[,1] == T),2])) > 0,x[which(x[,1] == T),2],0))
  names(bools) = c("False", "True")
  print(bools)
}
ClosePlotDevice = function(){
  while(!is.null(dev.list())){dev.off()}
}
UnifyData = function(df, vec, var.names = c("V1", "V2")) {
  for (i in 1:ncol(df)){
    if (i == 1){var1 = NULL
    var2 = NULL}
    
    var1 = c(var1, df[,i])
    var2 = c(var2, vec)
    
    if(i == ncol(df)){
      df.new = data.frame(V1 = var1, V2 = var2)
      names(df.new) = c(var.names)
      return(df.new)}
    
  }
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
#################################################################################################################
#                                              Importing data                                                   #
#################################################################################################################
setwd("F:/MG_Analysis_v2/Data/")
#setwd("C:/MG_Analysis_v2/Data/")
Meta = read.csv2("ANT28_MG_Meta.csv")

setwd("F:/MG_Analysis_v2/Data/AOM_Mapped_Datasets")
AOM = read.csv2(file = "AOM_mapped_counts_v3.csv", header = T, stringsAsFactors = F)
AOM[,18:39] = cpm(data = AOM[,18:39], gene_length = AOM$Gene_Length)

#################################################################################################################
#                                         CAZyme abundance vs. Temperature                                      #
#################################################################################################################
CAZy_main_counts = AOM[which(!is.na(AOM$CAZy_EC)),c(17:39)]
con = NULL
for (i in 2:ncol(CAZy_main_counts)){
  if (i == 1){con = NULL}  
  #print(names(CAZy_main_counts)[i])
  print(sum(CAZy_main_counts[,i]))
  con = c(con, sum(CAZy_main_counts[,i])) 
}

for (i in 1:(ncol(CAZy_main_counts)-1)){
  if (i == 1){

    df = CAZy_main_counts[,c(1,i+1)]
    df[,1] = as.factor(df[,1])
    names(df) = c("CAZy_EC", "Count")
    df.cazy.ec = summarise(group_by(df, CAZy_EC), sum = sum(Count))  
    names(df.cazy.ec) = c("CAZy_EC", colnames(CAZy_main_counts[i+1]))
  }else{
    df = CAZy_main_counts[,c(1,i+1)]
    df[,1] = as.factor(df[,1])
    names(df) = c("CAZy_EC", "Count")
    sum.cazy.ec = summarise(group_by(df, CAZy_EC), sum = sum(Count))  
    names(sum.cazy.ec) = c("CAZy_EC", colnames(CAZy_main_counts[i+1]))
    df.cazy.ec = merge(df.cazy.ec, sum.cazy.ec, by = "CAZy_EC", all = T)
  }
  if (i == (ncol(CAZy_main_counts)-1)){
    rownames(df.cazy.ec) = df.cazy.ec[,1]
    df.cazy.ec = df.cazy.ec[,-1]
  }
}
colSums(df.cazy.ec) == con
df.cazy.ec.nn = df.cazy.ec 
df.cazy.ec.t = t(df.cazy.ec)
df.cazy.ec.t = Normalize1(df.cazy.ec.t)
for (i in 1:ncol(df.cazy.ec.t)){
  if(i == 1) {
    vec.p = NULL
    vec.r = NULL
    vec.cazy = NULL
  }
  val = data.frame(Val = df.cazy.ec.t[,i], Temperature = Meta$Temperature_T)
  val = val[order(val$Temperature),] 
  names(val) = c(colnames(df.cazy.ec.t)[i], "Temperature")
  
  mod = lm(val[,1]~poly(val[,2],2,raw = T))
  s.mod = summary(mod)
  r = s.mod$adj.r.squared
  p = lmp(mod)
  
  vec.p = c(vec.p, p)
  vec.r = c(vec.r, r)
  vec.cazy = c(vec.cazy, colnames(val)[1])
  
  
if (p<=0.05 && r>=.3){
  par(mfrow = c(1,1))
  name = paste("CAZy_Temperature_", colnames(val)[1], ".png", sep = "")
  # Cairo(file = name,
  #       type = "png",
  #       width = 8,
  #       height = 6,
  #       bg = "white",
  #       units = "in",
  #       dpi = 96)
    fit.mod = predict(mod, interval = "conf")
  plot(val[,1]~ val[,2], type = "p", xlab = "Temperature", ylab = "Sequence abunadance [%]", pch = 21, bg = "black",
       main = paste(colnames(val)[1], "\n adj.r = ", round(r,3), ", pval = ", round(p,3), sep = ""))
  lines(fit.mod[,"fit"] ~ val[,2], col = "red", lty = 2, lwd = 2)
  lines(fit.mod[,"lwr"] ~ val[,2], col = "black", lty = 3, lwd = 1)
  lines(fit.mod[,"upr"] ~ val[,2], col = "black", lty = 3, lwd = 1)
#  dev.off()
}
  if (i == ncol(df.cazy.ec.t)){
    CAZy.mod.df = cbind.data.frame(vec.cazy, vec.r, vec.p)
    names(CAZy.mod.df) = c("CAZyID", "adj.r", "p.val")
  }
}

CAZy.mod.df$adj.p.val = p.adjust(CAZy.mod.df$p.val, method = "BH")
names(CAZy.mod.df) = c("CAZyID", "adj.r", "p.val", "adj.p.val")
table(CAZy.mod.df$adj.p.val <= 0.05)[2]/nrow(CAZy.mod.df)
setwd("F:/MG_Analysis_v2/")
write.csv2(file = "cazyMOD.csv", CAZy.mod.df)

#################################################################################################################
#                                             CAZyme abundance cluster                                          #
#################################################################################################################
summary(CAZy.mod.df$adj.p.val <= .05)
cazymes = c("AA", "CBM", "CE", "GH", "GT", "PL")

for(i in 1:length(cazymes)){
  if(i == 1){
    barplot.data = data.frame(matrix(nrow = 10, ncol = length(cazymes)))
    names(barplot.data) = cazymes
  }
  df = CAZy.mod.df[which(grepl(cazymes[i], CAZy.mod.df$CAZyID)),]
  table(df$adj.p.val <= 0.05) 
  barplot.data[,i] = c(SumBool(as.data.frame(table(abs(df$adj.r) > 0 & abs(df$adj.r) < .1 & df$adj.p.val <= .05)))[2],
                       SumBool(as.data.frame(table(abs(df$adj.r) > .1 & abs(df$adj.r) < .2 & df$adj.p.val <= .05)))[2],
                       SumBool(as.data.frame(table(abs(df$adj.r) > .2 & abs(df$adj.r) < .3 & df$adj.p.val <= .05)))[2],
                       SumBool(as.data.frame(table(abs(df$adj.r) > .3 & abs(df$adj.r) < .4 & df$adj.p.val <= .05)))[2],
                       SumBool(as.data.frame(table(abs(df$adj.r) > .4 & abs(df$adj.r) < .5 & df$adj.p.val <= .05)))[2],
                       SumBool(as.data.frame(table(abs(df$adj.r) > .5 & abs(df$adj.r) < .6 & df$adj.p.val <= .05)))[2],
                       SumBool(as.data.frame(table(abs(df$adj.r) > .6 & abs(df$adj.r) < .7 & df$adj.p.val <= .05)))[2],
                       SumBool(as.data.frame(table(abs(df$adj.r) > .7 & abs(df$adj.r) < .8 & df$adj.p.val <= .05)))[2],
                       SumBool(as.data.frame(table(abs(df$adj.r) > .8 & abs(df$adj.r) < .9 & df$adj.p.val <= .05)))[2],
                       SumBool(as.data.frame(table(abs(df$adj.r) > .9 & abs(df$adj.r) < 1 & df$adj.p.val <= .05)))[2])
}

df.cazy.norm = Normalize1(df.cazy.ec.t)
df.cazy.norm.sig = df.cazy.ec.t[,which(CAZy.mod.df$adj.p.val <= .05)]
mypal = c("#75a9b6","#dcba32","#d42c08")

cazy.dist = vegdist(t(df.cazy.norm.sig), method = "euc")
cazy.ward.cluster = hclust(cazy.dist, method = "ward.D2")
cut = cutree(cazy.ward.cluster, k = 3)

cazy.dist.c = vegdist(t(df.cazy.ec.t), method = "euc")
cazy.ward.cluster.c = hclust(cazy.dist.c, method = "ward.D2")
cut.2 = cutree(cazy.ward.cluster.c, k = 3)

setwd("F:/MG_Analysis_v2/Plots_v3/CAZymes")
summary(as.factor(cut))
par(mfrow=c(1,1))
Cairo(file = "CAZy_Cluster_tree.svg",
      type = "svg",
      width = 10,
      height = 10,
      bg = "white",
      units = "in",
      dpi = 96)
plot(as.phylo(cazy.ward.cluster),  type = "fan", cex = .6,  tip.color = mypal[cut], label.offset = .02, xpd = T)
dev.off()

setwd("F:/MG_Analysis_v2/Data/CAZy")
data = as.matrix(read.csv2("CAZy_temperature_lm.csv", row.names = 1))
bar.cols = colorRampPalette(brewer.pal(n = 11, name = "RdYlGn"))(10)
setwd("F:/MG_Analysis_v2/Diversity/NatCom")
d = 2
ClosePlotDevice()
Cairo(file = "Figure_4_Supp_Cluster_v3.svg",
      type = "svg",
      width = 10,
      height = 10,
      bg = "white",
      units = "in",
      dpi = 96)
par(mar=c(5.5,5.5,5.5,5.5))
par(xpd=T)
plot(as.phylo(cazy.ward.cluster),  type = "fan", cex = 1,  tip.color = mypal[cut], label.offset = .02)
par(xpd=F)
par(mar=c(5.1,5.1,4.1,2.1))
ClosePlotDevice()

ClosePlotDevice()
Cairo(file = "Figure_4_CAZy_Cluster_v4.pdf",
      type = "pdf",
      width = 12,
      height = 9,
      bg = "white",
      units = "in",
      dpi = 96)
layout(matrix(c(1,1,1,1,5,5,
                1,1,1,1,5,5,
                2,2,3,3,4,4,
                2,2,3,3,4,4), ncol = 6, nrow = 4, byrow = T))
par(mar=c(5.1,5.1,4.1,2.1))
barplot(as.matrix(barplot.data), space = c(1.5),las = 2, col = bar.cols, xlab = "Number of CAZyme families", main = "", xlim = c(0,50), hor = T, axes = F)
axis(1, label = seq(from = 0, to = 50, by = 10), at = seq(from = 0, to = 50, by = 10))

CAZy.Cluster.2 = df.cazy.norm.sig[,which(cut == 1)]
CAZy.Cluster.2.long = UnifyData(df = CAZy.Cluster.2, vec = Meta$Temperature_T)
CAZy.Cluster.2.long = CAZy.Cluster.2.long[order(CAZy.Cluster.2.long$V2),]

PlotLm(y = CAZy.Cluster.2.long$V1, x = CAZy.Cluster.2.long$V2,
       ylim = c(0,1), d = d, l1.col = "black", mod.pos = "front", lty = 1, predict = F,
       col = "#75a9b6aa", pch = 16, las = 1,
       ylab = "Normalized CAZyme abundance", xlab = "", plot.result = "top")

CAZy.Cluster.1 = as.data.frame(df.cazy.norm.sig[,which(cut == 2)])
CAZy.Cluster.1.long = UnifyData(CAZy.Cluster.1, Meta$Temperature_T)
CAZy.Cluster.1.long = CAZy.Cluster.1.long[order(CAZy.Cluster.1.long$V2),]
PlotLm(y = CAZy.Cluster.1.long$V1, x = CAZy.Cluster.1.long$V2,
       ylim = c(0,1), d = d, l1.col = "black", mod.pos = "front", lty = 1, predict = F,
       col = "#dcba32aa", pch = 16, las = 1,
       ylab = "", xlab = "Temperature [°C]", plot.result = "top")

CAZy.Cluster.3 = df.cazy.norm.sig[,which(cut == 3)]
CAZy.Cluster.3.long = UnifyData(CAZy.Cluster.3, Meta$Temperature_T)
CAZy.Cluster.3.long = CAZy.Cluster.3.long[order(CAZy.Cluster.3.long$V2),]
PlotLm(y = CAZy.Cluster.3.long$V1, x = CAZy.Cluster.3.long$V2,
       ylim = c(0,1), d = d, l1.col = "black", mod.pos = "front", lty = 1, predict = F,
       col = "#d42c08aa", pch = 16, las = 1,
       ylab = "", xlab = "", plot.result = "top")

plot(0,0, axes = F, xlab = "", ylab = "", type = "n") 
legend("center", legend = rev(c("0.0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1.0")), fill = rev(bar.cols), pt.cex = 1.3, cex = 1, bty = "n", title = expression(paste("r"^"2")))
dev.off()