#################################################################################################################
#                                     Cleaning Workspace & loading packages                                     #
#################################################################################################################
rm(list = ls(all=T))
library(oceanmap)#
library(ncdf4)#
library(raster)#
library(Cairo)#
library(cluster)#
library(ape)#

##### Readin meta data #####
setwd("F:/MG_Analysis_v2/Data")
Meta = read.csv2("ANT28_MG_Meta_v2.csv")
#################################################################################################################
#                                             Custom functions                                                  #
#################################################################################################################
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
ClosePlotDevice = function(){
  while(!is.null(dev.list())){dev.off()}
}

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
#                               Plotting Environmental data and station clusters                                #
#################################################################################################################
setwd("F:/MG_Analysis_v2/Plots/")
Cairo(file="Figure_1_T_Lat.svg", 
      type="svg",
      bg = "white",
      units="in", 
      width=6,
      height=4.5, 
      dpi=96)
par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))
plot(x = Meta$Latitude, y = Meta$Temperature_T, lty = 1,type = "l", col = "black", cex.lab = .7, cex.axis = .6, ylim = c(0,30),
     main = "", ylab = "Temperature [°C]", xlab = "Latitude [°N]", las = 1)
points(x = Meta$Latitude, y = Meta$Temperature_T, pch = 21,col = "black", bg = as.character(Meta$Prov_Col), cex = 1.2)
abline(v = c(-60.8, -57, -53, -43, -4.5, 7, 26, 38.5), col = "black", lty = 2)
ClosePlotDevice()

ANT28.Legend = rev(c("APLR", "ANTA", "SANT", "FKLD", "SATL", "WTRA", "NATR", "NAST", "NADR"))
Legend.colors = rev(c("#85b6bb", "#226874","#0a303b", "#78369f", "#d04131", "#d27120", "#f1b434", "#558a33", "#327250"))
Cairo(file="Figure_1_Part2_Province_legend_squares.svg", 
      type="svg",
      bg = "white",
      units="in", 
      width=5, 
      height=5, 
      dpi=96)
par(mfrow = c(1,1))
par(xpd=T)
plot(0,0, type = "n", axes = F, ylab = "", xlab = "")
legend("center", pch = 21, pt.bg = Legend.colors, legend = ANT28.Legend,  bty = "n", pt.cex = 1.5,cex = .9, horiz = F)
par(xpd=F)
ClosePlotDevice()

Cairo(file="Figure_Sx_Part1_T_Lat_Talk.svg", 
      type="svg",
      bg = "white",
      units="in", 
      width=12, 
      height=12, 
      dpi=96)
layout(matrix(c(1,1,2,2,9,9,
                3,3,4,4,9,9,
                5,5,6,6,9,9,
                7,7,8,8,10,10), nrow = 4, ncol = 6, byrow = T))
plot(x = Meta$Latitude, y = Meta$Temperature_T, lty = 2,type = "l", col = "black", cex.lab = .7, cex.axis = .6, ylim = c(0,30),
     main = "", ylab = "Temperature [°C]", xlab = "", las = 1)
points(x = Meta$Latitude, y = Meta$Temperature_T, pch = 21,col = "black", bg = as.character(Meta$Prov_Col), cex = 1.2)
tmp = data.frame(x = Meta$Latitude, y= Meta$Salinity_psu, col = Meta$Prov_Col)
tmp = tmp[which(complete.cases(tmp)),]
plot(x = tmp$x, y = tmp$y, lty = 2,type = "l", col = "black", cex.lab = .7, cex.axis = .6,
     main = "", ylab = "Salinity [psu]", xlab = "Latitude [°N]", las = 1)
points(x = tmp$x, y = tmp$y, pch = 21,col = "black", bg = as.character(tmp$col), cex = 1.6)

tmp = data.frame(x = Meta$Latitude, y= Meta$Phosphate_annual_mean_20m, col = Meta$Prov_Col)
tmp = tmp[which(complete.cases(tmp)),]
plot(x = tmp$x, y = tmp$y, lty = 2,type = "l", col = "black", cex.lab = .7, cex.axis = .6,
     main = "", ylab = "Annu. mean phosphate [µM]", xlab = "", las = 1)
points(x = tmp$x, y = tmp$y, pch = 21,col = "black", bg = as.character(tmp$col), cex = 1.6)

tmp = data.frame(x = Meta$Latitude, y = Meta$Nitrate_annual_mean_20m, col = Meta$Prov_Col)
tmp = tmp[which(complete.cases(tmp)),]
plot(x = tmp$x, y = tmp$y, lty = 2,type = "l", col = "black", cex.lab = .7, cex.axis = .6,
     main = "", ylab = "Annu. mean nitrate [µM]", xlab = "Latitude [°N]", las = 1)
points(x = tmp$x, y = tmp$y, pch = 21,col = "black", bg = as.character(tmp$col), cex = 1.6)

tmp = data.frame(x = Meta$Latitude, y= Meta$Biomass_production_ngC_L_h, col = Meta$Prov_Col)
tmp = tmp[which(complete.cases(tmp)),]
plot(x = tmp$x, y = tmp$y, lty = 2,type = "l", col = "black", cex.lab = .7, cex.axis = .6,
     main = "", ylab = "BPP [ng C L-1 h-1]", xlab = "Latitude [°N]", las = 1)
points(x = tmp$x, y = tmp$y, pch = 21, col = "black", bg = as.character(tmp$col), cex = 1.6)

tmp = data.frame(x = Meta$Latitude, y= Meta$POC, col = Meta$Prov_Col)
tmp = tmp[which(complete.cases(tmp)),]
plot(x = tmp$x, y = tmp$y, lty = 2,type = "l", col = "black", cex.lab = .7, cex.axis = .6,
     main = "", ylab = "POC [mg/L]", xlab = "Latitude [°N]", las = 1)
points(x = tmp$x, y = tmp$y, pch = 21,col = "black", bg = as.character(tmp$col), cex = 1.6)

tmp = data.frame(x = Meta$Latitude, y= Meta$Chla_ugL, col = Meta$Prov_Col)
tmp = tmp[which(complete.cases(tmp)),]
plot(x = tmp$x, y = tmp$y, lty = 2,type = "l", col = "black", cex.lab = .7, cex.axis = .6,
     main = "", ylab = "Chl a [µg/L]", xlab = "Latitude [°N]", las = 1)
points(x = tmp$x, y = tmp$y, pch = 21,col = "black", bg = as.character(tmp$col), cex = 1.6)
tmp = data.frame(x = Meta$Latitude, y= Meta$Bateria_ml, col = Meta$Prov_Col)
tmp = tmp[which(complete.cases(tmp)),]
tmp$y = tmp$y/1000000
plot(x = tmp$x, y = tmp$y, lty = 2,type = "l", col = "black", cex.lab = .7, cex.axis = .6,
     main = "", ylab = "Bacterial cells x 107 [ml-1]", xlab = "Latitude [°N]", las = 1)
points(x = tmp$x, y = tmp$y, pch = 21,col = "black", bg = as.character(tmp$col), cex = 1.6)
meta.daisy = Meta[,c(12,13,19)] #16, bac;  | 19 chl | 26 POC |
meta.daisy2 = meta.daisy
for (i in 1:ncol(meta.daisy2)){
  t = meta.daisy2[,i]
  t[is.na(t)] = mean(meta.daisy2[,i], na.rm =T)
  meta.daisy2[,i] = t
}
meta.daisy2 = Normalize1(meta.daisy2)
meta.daisy2[is.na(meta.daisy)] = NA
rownames(meta.daisy2) = Meta$Station
clust = hclust(daisy(meta.daisy2, metric = "euc"), method = "ward.D2")
plot(as.phylo(clust), type = "phy", xpd = T)
tiplabels(pch = 21,
          bg = as.character(Meta$Prov_Col), cex = 1.5, offset = .22, xpd = T)
#par(xpd=T)
plot(0,0, type = "n", axes = F, ylab = "", xlab = "")
legend("center", pch = 21, pt.bg = Legend.colors, legend = ANT28.Legend,  bty = "n", pt.cex = 1.5,cex = .9, horiz = F)
#par(xpd=F)
ClosePlotDevice()

##### Plotting Province Legend #####
ANT28.Legend = rev(c("APLR", "ANTA", "SANT", "FKLD", "SATL", "WTRA", "NATR", "NAST", "NADR"))
Legend.colors = rev(c("#85b6bb", "#226874","#0a303b", "#78369f", "#d04131", "#d27120", "#f1b434", "#558a33", "#327250"))
Cairo(file="Figure_1_Part2_Province_legend_squares.svg", 
      type="svg",
      bg = "white",
      units="in", 
      width=5, 
      height=5, 
      dpi=96)
par(mfrow = c(1,1))
par(xpd=T)
plot(0,0, type = "n", axes = F, ylab = "", xlab = "")
legend("center", pch = 21, pt.bg = Legend.colors, legend = ANT28.Legend,  bty = "n", pt.cex = 1.5,cex = .9, horiz = F)
par(xpd=F)
dev.off()

#### Download AQUA MODIS Chlorophyll data from 2012 from: 
#### https://oceandata.sci.gsfc.nasa.gov/directaccess/MODIS-Aqua/Mapped/Annual/4km/chlor_a/
#### Regestration and login required!
setwd("F:/MG_Analysis_v2/Data/MODIS_CHLa_data")
chl = nc_open("A20120012012366.L3m_YR_CHL_chlor_a_4km.nc")
chl.dat.raster = nc2raster(chl, "chlor_a", lonname="lon", latname="lat", date=T)

sst.flip = flip(sst.dat.raster, "y")
chl.flip = flip(chl.dat.raster, "y")
lon = c(-75,5)
lat = c(-75,60)
### Transect Map
crop.vals = c(lon, lat)
sst.crop = raster::crop(sst.flip, extent(crop.vals)) 
chl.crop = raster::crop(chl.flip, extent(crop.vals)) 

##### Plotting Oceanmap with Chl a backgroung #####
##### Save Plot from preview windos as svg!!! #####
chl.pal = colorRampPalette(c("white", "darkgreen"))(300)

par(mfrow=c(1,1))
dev.off()
v(chl.crop, cbpos = "b", pal = chl.pal, zlim = c(0,1.5), cb.xlab = expression("Annual chlorophyll-a (mg m"^-3*")"),
  bwd = 0.01, grid = F, replace.na = F, col.bg = "white", border = "#504f4f")
abline(v = c(-50, -25, 0), col = "black", lty = c(2,2,1))
abline(h = c(-50, -25, 0, 25, 50), col = "black", lty = c(2,2,1,2,2))
points(y=Meta$Latitude, x = Meta$Longitude, pch = 21, col = "black", bg = as.character(Meta$Prov_Col), cex = 1.9)
text(y=Meta$Latitude, x = Meta$Longitude, Meta$Station, cex=.45, col = "white")
text(y=Meta$Latitude, x = Meta$Longitude, Meta$Station, col = "white", cex =.45)