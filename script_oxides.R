###Script for oxides estimation from spectra###

#### Install and load packages
if(!require(prospectr)) install.packages(prospectr) 
if(!require(corrplot)) install.packages(corrplot) 
if(!require(caret)) install.packages(caret)
if(!require(raster)) install.packages(raster)
library(raster)
require(raster)

## Set working directory
setwd("C:/Users/neliq/Google Drive/dados/data")
list.files()

## Load the data and the satellite image (SYSI)
data <- read.csv("dados_spectra.csv", h=TRUE, sep=";")

rownames(data) <- data[,1]  ##Convert the first column into rownames

spectra <- data[,-1:-3]
colnames(spectra) <- seq(from = 350, to = 2500, by=1)

SYSI <- stack("SYSI.tif")   #Synthetic Soil Image with six spectral bands
plot(SYSI[[1]])
SYSI[SYSI == 0] <- NA

## Pre-processing by Savitzky and Golay (1968)
## Plot raw spectra

matplot(as.numeric(colnames(spectra)), t(spectra[1:4,]), type = "l",
        xlab = "wavelenght (nm)",
        ylab = "Reflectance factor")

### Smooth spectrum without deriving

spectra.process1 <- savitzkyGolay(spectra, m = 0, p = 3, w = 35)  #moving windows, it must be odd
matplot(as.numeric(colnames(spectra.process1)), t(spectra.process1[1:4,]), type = "l",
        xlab = "wavelenght (nm)",
        ylab = "Reflectance factor")

### Kubelka-Munk transformation and second derivative

KM = ((1-spectra.process1)^2)/(2*spectra.process1) # from smooth spectra
matplot(as.numeric(colnames(KM)), t(KM[1:5,]), type = "l",
        xlab = "wavelenght (nm)",
        ylab = " ")

sec.deriv = savitzkyGolay(KM, m = 2, p = 3, w = 25)
matplot(as.numeric(colnames(sec.deriv[,100:250])), t(sec.deriv[1:4,100:250]), type = "l",
        xlab = "wavelenght (nm)",
        ylab = " ")

## Estimates of AHm and AGm

### AHm
Hm.max = apply(data.frame(sec.deriv[ ,157:202]), MARGIN = 1, FUN = max)  #The sec.deriv[,157:202] 
#corresponds to the amplitude between 415 and 445 nm

Hm.min = apply(data.frame(sec.deriv[ ,157:202]), MARGIN = 1, FUN = min) 

AHm = as.data.frame((Hm.max^2)^(1/2))+((Hm.min^2)^(1/2))
AHm = as.data.frame(Hm.max-Hm.min)
colnames(AHm) <- "AHm"
head(AHm)

### AGt
Gt.max = apply(data.frame(sec.deriv[ ,37:77]), MARGIN = 1, FUN = max) #The sec.deriv[,37:67] 
#corresponds to the amplitude between 415 and 445 nm

Gt.min = apply(data.frame(sec.deriv[ ,37:67]), MARGIN = 1, FUN = min)

AGt = as.data.frame((Gt.max^2)^(1/2))+((Gt.min^2)^(1/2))
AGt1 = as.data.frame(Gt.max-Gt.min)
colnames(AGt) <- "AGt"
head(AGt)

Hm.Gt <- cbind(AHm, AGt)
head(Hm.Gt)

ratio <- as.data.frame(Hm.Gt$AHm / (Hm.Gt$AHm + Hm.Gt$AGt))
colnames(ratio) <- "ratio"
rownames(ratio) <- rownames(AHm)
head(ratio)

### Plot the SYSI and the sampling points
dat2 <- cbind(data[,2:3], ratio)
coordinates(dat2) <- ~X+Y
plotRGB(SYSI, r=3, g=2, b=1, stretch = "lin")
plot(ratio, add=T, col = "red", pch = 20)
spplot(dat2)

### Hm and Gt according to Fernandes et al. (2004)

Hm_gkg = -1.6 + (36320*AHm)  #(R2 = 0.94)
Gt_gkg = 5.7 + (18607*AGt)   #(R2 = 0.63)


dat_Hm_Gt <- as.data.frame(cbind(data[,2:3], Hm_gkg, Gt_gkg))
colnames(dat_Hm_Gt) <- c("X", "Y", "Hm_g_kg", "Gt_g_kg")
coordinates(dat_Hm_Gt)<- ~X+Y
spplot(dat_Hm_Gt, colorkey = TRUE)


### Hematite and goethite ratio (Hm/(Hm+Gt)) proposed by Fernandes et al. (2004)
Hm_Gt_calculated = -0.059 + (1.506*AHm)/(AHm+AGt)
Hm.Gt.ratio <- cbind(data[,2:3], Hm_Gt_calculated)
colnames(Hm.Gt.ratio) <- c("X", "Y", "ratioHm_Gt")
coordinates(Hm.Gt.ratio) <- ~X+Y
spplot(Hm.Gt.ratio, colorkey=TRUE)

## Niodi Index proposed by Viscarra Rossel et al. (2010)

lenght = 367:2483
CR = continuumRemoval(spectra.process1, wav = lenght)
matplot(as.numeric(colnames(CR)), t(CR[1:5,]), type = "l",
        xlab = "wavelenght (nm)",
        ylab = "Continuum Removal")

D = 1 - CR
D880 = apply(data.frame(D[ ,499:569]), MARGIN = 1, FUN = max)
D920 = apply(data.frame(D[ ,549:609]), MARGIN = 1, FUN = max)

NIODI = as.data.frame((D920 - D880)/(D880 + D920))
colnames(NIODI) <- "NIODI"


dat3 <- cbind(data[,2:3], NIODI)
coordinates(dat3) <- ~X+Y
spplot(dat3, colorkey = TRUE)



## Kaolinite and Gibbsite from Contiuum removal 

### Calculate minumum and maximum reflectance between  and 2293 nm

CR_min = apply(data.frame(D[ ,1835:1898]), MARGIN = 1, FUN = min)
KT_max = apply(data.frame(D[ ,1766:1858]), MARGIN = 1, FUN = max)
GB_max = apply(data.frame(D[ ,1858:1927]), MARGIN = 1, FUN = max)


Kt_index = as.data.frame(KT_max - CR_min)
Gb_index = as.data.frame(GB_max - CR_min)

Kt_Gb <- cbind(data[,2:3], Kt_index, Gb_index)
colnames(Kt_Gb) <- c("X", "Y", "Kt", "Gb")
coordinates(Kt_Gb) <- ~X+Y
spplot(Kt_Gb, colorkey = TRUE)

#Kt/(Kt+Gb)
ratioKt_Gb <- (Kt_Gb$Kt/(Kt_Gb$Kt+Kt_Gb$Gb))
ratioKt_Gb_withXY <- cbind(data[,2:3],ratioKt_Gb)
colnames(ratioKt_Gb_withXY) <- c("X", "Y", "ratioKt_Gb")
coordinates(ratioKt_Gb_withXY) <- ~X+Y
spplot(ratioKt_Gb_withXY, colorkey = TRUE)

## Join all estimated values together

data.mineralogy <- cbind(data[,2:3], ratio, NIODI, Hm_gkg, Gt_gkg, Hm_Gt_calculated, ratioKt_Gb)
colnames(data.mineralogy) <- c("X", "Y", "ratioHm_Gt", "NIODI", "Hm_g_kg", "Gt_g_kg",
                               "ratioHm_Gt_calculated", "ratioKt_Gb")
head(data.mineralogy)

## Spatial prediction of the estimates mineralogical indexes and contents
coordinates(data.mineralogy) <- ~X+Y  #make spatial object

SYSI.df <- extract(SYSI, data.mineralogy)  #estract SYSI values for each soil sample
head(SYSI.df)

dat1 <- cbind(as.data.frame(data.mineralogy), SYSI.df) #join data and SYSI values in a single table
head(dat1)

### Pearson?s correlation between mineralogical index and SYSI bands
corrplot.mixed(cor(dat1[,-1:-2]), number.cex = 0.7, #correlation between SYSI 
               tl.cex = 0.8)  #and mineralogical values



## Mapping mineralogical indexes and amounts using satellite images

### Construct the regression three model via Cubist for Hm/(Hm+Gt) ratio
ratioHm.Gt <- train(dat1[,9:14], dat1$ratioHm_Gt, method = "cubist", 
                    trControl = trainControl(method = "cv", number = 10))
ratioHm.Gt

### Predict on the satellite 
ratioHm.Gt.map <- raster::predict(SYSI, ratioHm.Gt)
plot(ratioHm.Gt.map)

