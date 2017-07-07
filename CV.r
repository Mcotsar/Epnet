###### HOW TO CALCULATE THE CV IN AMPHORAE PRODUCTION ####

library(plyr) #works to count

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

#library raster 

library(raster)

mydata <- function(mean, sd){
      (sd/mean)*100
}

mydata(mean = , sd = )
































