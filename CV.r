###### HOW TO CALCULATE THE CV IN AMPHORAE PRODUCTION ####

library(plyr) #works to count

library(raster)#calculate the coefficient of variation 


#from the website https://alstatr.blogspot.com.es/2013/06/measure-of-relative-variability.html

CV <- read.csv('drespaper.csv', header=T, sep=",")

##if you want to choose the type of Dressel to do the analyse. We will choose three types of Dressel 20 amphorae
CV = subset(CV, type %in% c("Dressel C","Dressel D","Dressel E"))
CV <- CV[,5:13]
#to count the data
count(CV)

##we have to calculate the mean from each workshop with the measurement protruding_rim as example

#Malpica

MN1 = subset(CV, site=="malpica")
mpl = (MN1$protruding_rim)
mean(mpl)

MN1 = subset(CV, site=="malpica")
mpl = (MN1$exterior_diam)
mean(mpl)

#Delicias

MN2= subset(CV, site=="delicias")
dlc = (MN2$protruding_rim)
mean(dlc)

#Parlamento

MN3= subset(CV, site=="parlamento")
prl = (MN3$protruding_rim)
mean(prl)

#Cerro del Belén

MN4= subset(CV, site=="belen")
bln = (MN4$protruding_rim)
mean(bln)

#Villaseca

MN5= subset(CV, site=="villaseca")
bln = (MN5$protruding_rim)
mean(bln)


##Calculate the standard deviation from each workshop with the measurement protruding_rim as example

#Malpica

SD1 = subset(CV, site=="malpica")
mpl1 = (SD1$protruding_rim)
sd(mpl1)

SD1 = subset(CV, site=="malpica")
mpl1 = (SD1$exterior_diam)
sd(mpl1)

#Delicias

SD2= subset(CV, site=="delicias")
dlc1 = (SD2$protruding_rim)
sd(dlc1)

#Parlamento

SD3= subset(CV, site=="parlamento")
prl1 = (SD3$protruding_rim)
sd(prl1)

#Cerro del Belén

SD4= subset(CV, site=="belen")
bln1 = (SD4$protruding_rim)
sd(bln1)

#Villaseca

SD5= subset(CV, site=="villaseca")
bln1 = (SD5$protruding_rim)
sd(bln1)

##Calculate the coefficient of variation with the measurement protruding_rim as example


mydata <- function(mean, sd){
      (sd/mean)*100
}

mydata(sd = , mean = )


###EXAMPLE WITH MALPICA AND EXTERIOR_DIAM (MORE STANDARIZED MEASUREMENT)

MN1 = subset(CV, site=="malpica")
mpl = (MN1$exterior_diam)
mean(mpl)

SD1 = subset(CV, site=="malpica")
mpl1 = (SD1$exterior_diam)
sd(mpl1)

mydata <- function(mean, sd){
      (sd/mean)*100
}

mydata(sd = , mean = )


#SOME RESULTS WITH EXTERIOR_DIAM MEASUREMENT###

#MALPICA = sd(9.83), mean (166.05) -------> CV = 5.9199

#DELICIAS = sd(8.52), mean (172.084) -------> CV = 4.95407

#PARLAMENTO = sd(11.64), mean (163.8095) -------> CV = 7.106227

#BELEN = sd(12.37), mean (171.1818) -------> CV = 7.226235

#VILLASECA = sd(12.24), mean (160.2075) -------> CV = 7.64

























