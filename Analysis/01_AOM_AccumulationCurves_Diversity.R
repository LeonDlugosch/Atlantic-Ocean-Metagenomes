#################################################################################################################
#                                     Cleaning Workspace & loading packages                                     #
#################################################################################################################
rm(list = ls(all=T))
library(Cairo)
library(dplyr)
library(drc)
library(picante)
library(rtk)
library(iNEXT)
library(RAM)
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

#################################################################################################################
#                                               Importing data                                                  #
#################################################################################################################
setwd("F:/MG_Analysis_v2/Data/")
Meta = read.csv2("ANT28_MG_Meta_v2.csv")
lm.cols = data.frame(Temp = Meta$Temperature_T, Col = Meta$Color)
##### Reading AOM #####
setwd("F:/MG_Analysis_v2/Data/AOM_Mapped_Datasets")
data = read.csv2("AOM_mapped_counts_v3.csv")
data = data[,-2]
data.cpm = read.csv2("AOM_mapped_counts_CPM_v3.csv")
data.cpm.red = data.cpm[which(!is.na(data.cpm$Domain) & !is.na(data.cpm$Knr)),] 
all(data$SeqID == data.cpm$SeqID)

#################################################################################################################
#                                               Generating datasets                                             #
#################################################################################################################

data.red = data[which(!is.na(data$Domain) & !is.na(data$Knr)),] 
species = SummarizeDataset(data.red[,17:38], by = data.red$Species)
knr = SummarizeDataset(data.red[,17:38], by = data.red$Knr)

#################################################################################################################
#                           Boootstrapped rarefication and Diversity analysis                                   #
#################################################################################################################
#### This may take some time
plot.RTK = T
boots = 99
depth = 2000000

for(k in 1){
  for (i in 1:boots){
    ##### nr dataset ##### 
    if(i == 1){
      seq.rich = NULL
      seq.shan = NULL
      seq.shan.en = NULL
      seq.simp = NULL
      seq.simp.en = NULL
      seq.chao = NULL
      seq.even = NULL
    } 
    seq.rtk = rtk(seq.red, depth = depth, ReturnMatrix = 1)
    seq = as.data.frame(seq.rtk$raremat[[1]])
    seq.rich = cbind(seq.rich, specnumber(t(seq))) 
    seq.shan = cbind(seq.shan, diversity(t(seq)))
    seq.shan.en = cbind(seq.shan.en, exp(seq.shan))
    seq.simp = cbind(seq.simp, diversity(t(seq), index = "simpson"))
    seq.simp.en = cbind(seq.simp.en, diversity(t(seq), index = "invsimpson"))
    seq.chao = cbind(seq.chao, estimateR(x = t(seq))[2,])
    seq.even = cbind(seq.even, seq.shan/log(seq.rich))
    
    
    if(i == boots){
      seq.diversity = data.frame(SeqRichness = rowMeans(seq.rich),
                                 SeqShannon = rowMeans(seq.shan),
                                 SeqShannon.EN = rowMeans(seq.shan.en),
                                 SeqSimpson = rowMeans(seq.simp),
                                 SeqSimpson.EN = rowMeans(seq.simp.en),
                                 SeqChaoI = rowMeans(seq.chao),
                                 SeqEveness = rowMeans(seq.even))
    }
  }
  ##### taxonomic dataset ##### 
  for (i in 1:boots){
    if(i == 1){
      spe.rich = NULL
      spe.shan = NULL
      spe.shan.en = NULL
      spe.simp = NULL
      spe.simp.en = NULL
      spe.chao = NULL
      spe.even = NULL
    }  
    
    species.rtk = rtk(species.red, depth = depth, ReturnMatrix = 1)
    spe = as.data.frame(species.rtk$raremat[[1]])
    
    spe.rich = cbind(spe.rich, specnumber(t(spe))) 
    spe.shan = cbind(spe.shan, diversity(t(spe)))
    spe.shan.en = cbind(spe.shan.en, exp(spe.shan))
    spe.simp = cbind(spe.simp, diversity(t(spe), index = "simpson"))
    spe.simp.en = cbind(spe.simp.en, diversity(t(spe), index = "invsimpson"))
    spe.chao = cbind(spe.chao, estimateR(x = t(spe))[2,])
    spe.even = cbind(spe.even, spe.shan/log(spe.rich))
    if(i == boots){
      spe.diversity = data.frame(SeqRichness = rowMeans(spe.rich),
                                 SeqShannon = rowMeans(spe.shan),
                                 SeqShannon.EN = rowMeans(spe.shan.en),
                                 SeqSimpson = rowMeans(spe.simp),
                                 SeqSimpson.EN = rowMeans(spe.simp.en),
                                 SeqChaoI = rowMeans(spe.chao),
                                 SeqEveness = rowMeans(spe.even))
    }
  }
  ##### KO dataset ##### 
  for (i in 1:boots){
    if(i == 1){
      knr.rich = NULL
      knr.shan = NULL
      knr.shan.en = NULL
      knr.simp = NULL
      knr.simp.en = NULL
      knr.chao = NULL
      knr.even = NULL
    }  
    knr.rtk = rtk(knr.red, depth = 2000000, ReturnMatrix = 1)
    knr = as.data.frame(knr.rtk$raremat[[1]])
    
    knr.rich = cbind(knr.rich, specnumber(t(knr))) 
    knr.shan = cbind(knr.shan, diversity(t(knr)))
    knr.shan.en = cbind(knr.shan.en, exp(knr.SID))
    knr.simp = cbind(knr.simp, diversity(t(knr), index = "simpson"))
    knr.simp.en = cbind(knr.simp.en, diversity(t(knr), index = "invsimpson"))
    knr.chao = cbind(knr.chao, estimateR(x = t(knr))[2,])
    knr.even = cbind(knr.even, knr.shan/log(knr.rich))
    if(i == boots){
      knr.diversity = data.frame(KnrRichness = rowMeans(knr.rich),
                                 KnrShannon = rowMeans(knr.shan),
                                 KnrShannon.EN = rowMeans(knr.shan.en),
                                 KnrSimpson = rowMeans(knr.simp),
                                 KnrSimpson.EN = rowMeans(knr.simp.en),
                                 KnrChaoI = rowMeans(knr.chao),
                                 KnrEveness = rowMeans(knr.even))
    }
  }
}
##### Writing Diversity Data to .csv #####
setwd("F:/MG_Analysis_v2/Data/Supplement/")
Diversity = data.frame(Species_Richness = spe.SPR, Species_TSH = spe.TSH, Species_Evenness = spe.EVE,
                       KOG_Richness = knr.SPR, KOG_TSH = knr.TSH, KO_Evenness = knr.EVE,
                       nrAOM_seq_Richness = seq.SPR, nrAOM_seq_TSH = seq.TSH, nrAOM_seq_Evenness = seq.EVE)

write.csv2(Diversity, "TableS3_Diversity.csv")

#################################################################################################################
#                                             Collectors-curve plots                                            #
#################################################################################################################
#### This may take some time

setwd("F:/MG_Analysis_v2/Plots_v3")
if (plot.RTK == T){
  Cairo(file="ANT28_Collectors_curves.svg", 
        type="svg",
        bg = "white",
        units="in", 
        width=8, 
        height=12, 
        dpi=96)
  par(mfrow = c(3,1), cex.axis = 1, cex.lab = 1)
  rtk::collectors.curve(seq.rtk, main = "", times = 100, bin = 1,
                        xlab ="", ylab = expression(paste("nr-sequence richness (x10"^"6",")")), col = "grey", col2 = "black",
                        axes = F) 
  axis(1)
  axis(2, labels = c("0.0", "0.5", "1.0", "1.5", "2.0", "2.5"), at = c(0, 0.5e6, 1e6, 1.5e6, 2e6, 2.5e6), las = 2)
  box()
  rtk::collectors.curve(knr.rtk, main = "", times = 100, bin = 1,
                        xlab ="", ylab = expression(paste("KO richness (x10"^"3",")")), col = "grey", col2 = "black",
                        axes = F) 
  axis(1)
  axis(2, labels = c("0.0", "3.0", "6.0", "9.0", "12.0"), at = c(0, 3000, 6000, 9000, 12000), las = 2)
  box()
  rtk::collectors.curve(species.rtk, main = "", times = 100, bin = 1,
                        xlab ="No. of samples", ylab = expression(paste("Taxon richness (x10"^"4",")")), col = "grey", col2 = "black",
                        axes = F) 
  axis(1)
  axis(2, labels = c("0.0", "0.5", "1.0", "1.5", "2.0"), at = c(0, 5000, 10000, 15000, 20000), las = 2)
  box()
  dev.off()
}