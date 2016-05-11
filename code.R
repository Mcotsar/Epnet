library(ggplot2)
library(gridExtra)    
library(plyr)
library(MASS)
library(caret)

#####################CHECK OLD RESULTS##################################################

myData <- read.csv('all.csv', header=T, sep=";")
myData <- subset(myData, TYPE=="2")
myData <- myData[,5:13]

##ready for PCA!!
logData <- log(myData[,1:8])
pcaResults <- princomp(logData, center=T, scale=T)
plot(pcaResults)
pcaValues <- as.data.frame(pcaResults$scores)
pcaValues$site <- myData$site
ggplot(pcaValues, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point()
ggplot(pcaValues, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point() + facet_wrap(~site,ncol=1)

#

sampleSize = min(count(myData,'site')$freq)

amph1 <- subset(myData, site=="delicias")
sample1 <- amph1[sample(nrow(amph1), sampleSize),]
amph2 <- subset(myData, site=="malpica")
sample2 <- amph2[sample(nrow(amph2), sampleSize),]
amph3 <- subset(myData, site=="belen")
sample3 <- amph3[sample(nrow(amph3), sampleSize),]

sample <- rbind(sample1,sample2)
sample <- rbind(sample, sample3)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- sample$site

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1,1)/3)
predQda <- predict(amphDA, pcaValuesSample)
sample$probDelicias<- predQda$posterior[,"Delicias"]
sample$probDelicias<- predQda$posterior[,"Las Delicias"]
sample$probMalpica <- predQda$posterior[,"Malpica"]
sample$probBelen <- predQda$posterior[,"belén"]
pcaValuesSample$class <- predQda$class

#CONFUSION

confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)

svg('puntitos.svg')    
ggplot(pcaValuesSample, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point() + facet_wrap(~correcto, ncol=1) + ylim(c(-0.3,0.3))
dev.off()

foo <- subset(pcaValuesSample, Comp.1>-1.2)
foo <- subset(foo, Comp.2>-2.1)

svg('fig_dist.svg')    
ggplot(foo, aes(x=Comp.1, y=Comp.2)) + geom_density2d(aes(col=site), alpha=0.3) + geom_point(aes(col=site), size=2) + facet_wrap(~site, ncol=1) + theme(legend.position='none')
dev.off()

####### peer to peer #############

## init matrix
distMetrics <- matrix(0, nrow=3, ncol=3)
rownames(distMetrics) <- c('belén','Las Delicias','Malpica')
colnames(distMetrics) <- c('belén','Las Delicias','Malpica')

    
belen <- subset(myData, site=='belén')    
delicias <- subset(myData, site=='Las Delicias')    
malpica <- subset(myData, site=='Malpica')    

# 1 - Malpica-belen

pairData <- rbind(malpica,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- malpica[sample(nrow(malpica), sampleSize),]
sample2 <- belen[sample(nrow(belen), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['Malpica','belén'] <- conf$overall['Accuracy']
distMetrics['belén','Malpica'] <- conf$overall['Accuracy']


# 2 - Malpica-Las Delicias

pairData <- rbind(malpica,delicias)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- delicias[sample(nrow(delicias), sampleSize),]
sample2 <- malpica[sample(nrow(malpica), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['Malpica','Las Delicias'] <- conf$overall['Accuracy']
distMetrics['Las Delicias','Malpica'] <- conf$overall['Accuracy']


# 3 - Las Delicias-belén

pairData <- rbind(delicias,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- delicias[sample(nrow(delicias), sampleSize),]
sample2 <- belen[sample(nrow(belen),sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['Las Delicias','belén'] <- conf$overall['Accuracy']
distMetrics['belén','Las Delicias'] <- conf$overall['Accuracy']

distMetrics


#############################CHECK WITH NEW RESULTS############################################


myData <- read.csv('dres.csv', header=T, sep=",")
myData= subset(myData, type %in% c("Dressel C","Dressel D","Dressel E","Dressel G","Dressel H"))
myData <- myData[,5:13]
myData=myData[myData$site=="delicias" | myData$site=="malpica" | myData$site== "belen",]  
myData$site = factor(myData$site)

##ready for PCA!!
logData <- log(myData[,1:8])
pcaResults <- princomp(logData, center=T, scale=T)
plot(pcaResults)

pcaValues <- as.data.frame(pcaResults$scores)
pcaValues$site <- myData$site
ggplot(pcaValues, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point()
ggplot(pcaValues, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point() + facet_wrap(~site,ncol=1)


###discriminant###########

sampleSize = min(count(myData,'site')$freq)

amph1 <- subset(myData, site=="delicias")
sample1 <- amph1[sample(nrow(amph1), sampleSize),]
amph2 <- subset(myData, site=="malpica")
sample2 <- amph2[sample(nrow(amph2), sampleSize),]
amph3 <- subset(myData, site=="belen")
sample3 <- amph3[sample(nrow(amph3), sampleSize),]

sample <- rbind(sample1,sample2)
sample <- rbind(sample, sample3)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- sample$site

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1,1)/3)
predQda <- predict(amphDA, pcaValuesSample)
sample$probDelicias<- predQda$posterior[,"delicias"]
sample$probMalpica <- predQda$posterior[,"malpica"]
sample$probBelen <- predQda$posterior[,"belen"]
pcaValuesSample$class <- predQda$class

###CONFUSION###


confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)

svg('puntitos2.svg')    
ggplot(pcaValuesSample, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point() + facet_wrap(~correcto, ncol=1) + ylim(c(-0.3,0.3))
dev.off()

foo <- subset(pcaValuesSample, Comp.1>-1.2)
foo <- subset(foo, Comp.2>-2.1)

svg('fig_dist2.svg')    
ggplot(foo, aes(x=Comp.1, y=Comp.2)) + geom_density2d(aes(col=site), alpha=0.3) + geom_point(aes(col=site), size=2) + facet_wrap(~site, ncol=1) + theme(legend.position='none')
dev.off()


####### peer to peer #############

## init matrix
distMetrics <- matrix(0, nrow=3, ncol=3)
rownames(distMetrics) <- c('belen','delicias','malpica')
colnames(distMetrics) <- c('belen','delicias','malpica')

belen <- subset(myData, site=='belen')    
delicias <- subset(myData, site=='delicias')    
malpica <- subset(myData, site=='malpica')    

# 1 - Malpica-belen

pairData <- rbind(malpica,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- malpica[sample(nrow(malpica), sampleSize),]
sample2 <- belen[sample(nrow(belen), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['malpica','belen'] <- conf$overall['Accuracy']
distMetrics['belen','malpica'] <- conf$overall['Accuracy']


# 2 - Malpica-Las Delicias

pairData <- rbind(malpica,delicias)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- delicias[sample(nrow(delicias), sampleSize),]
sample2 <- malpica[sample(nrow(malpica), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['malpica','delicias'] <- conf$overall['Accuracy']
distMetrics['delicias','malpica'] <- conf$overall['Accuracy']


# 3 - Las Delicias-belén

pairData <- rbind(delicias,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- delicias[sample(nrow(delicias), sampleSize),]
sample2 <- belen[sample(nrow(belen),sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['delicias','belen'] <- conf$overall['Accuracy']
distMetrics['belen','delicias'] <- conf$overall['Accuracy']

distMetrics

################DOING THE SAME BUT WITH ANOTHER WORKSHOP:PARLAMENTO###############################################################################################################

myData <- read.csv('dres.csv', header=T, sep=",")
myData= subset(myData, type %in% c("Dressel C","Dressel D","Dressel E","Dressel G","Dressel H"))
myData <- myData[,5:13]


sampleSize = min(count(myData,'site')$freq)

amph1 <- subset(myData, site=="delicias")
sample1 <- amph1[sample(nrow(amph1), sampleSize),]
amph2 <- subset(myData, site=="malpica")
sample2 <- amph2[sample(nrow(amph2), sampleSize),]
amph3 <- subset(myData, site=="belen")
sample3 <- amph3[sample(nrow(amph3), sampleSize),]
amph4 <- subset(myData, site=="parlamento")
sample4 <- amph4[sample(nrow(amph4), sampleSize),]

sample <- rbind(sample1,sample2)
sample <- rbind(sample, sample3)
sample <- rbind(sample, sample4)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)

# plot to check the relevance of the first 2 PC'S
plot(pcaResultsSample)


# get the scores of the data
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
# put type in pcaValues
pcaValuesSample$site <- sample$site
ggplot(pcaValuesSample, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point()
ggplot(pcaValuesSample, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point() + facet_wrap(~site,ncol=1)


amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1,1,1)/4)
predQda <- predict(amphDA, pcaValuesSample)
    
sample$probDelicias<- predQda$posterior[,"delicias"]
sample$probMalpica <- predQda$posterior[,"malpica"]
sample$probBelen <- predQda$posterior[,"belen"]
sample$probParlamento <- predQda$posterior[,"parlamento"]
pcaValuesSample$class <- predQda$class

#g1 <- ggplot(pcaValuesSample, aes(x=Comp.1, y=Comp.2, colour=factor(class))) + geom_point() + facet_grid(~site) + ggtitle("pca1_2")

#ggplot(pcaValuesSample, aes(x=Comp.1, y=Comp.2, col=interaction(site,class), label=site)) + geom_text(size=5) + theme_bw() + theme(legend.position="top")

confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)

correct <- subset(pcaValuesSample, site==class)
incorrect <- subset(pcaValuesSample, site!=class)
correct$predict <- "pred. correcta"
incorrect$predict <- "pred. incorrecta"

result <- rbind(correct, incorrect)
result <- subset(result, Comp.1>-1.2)
result <- subset(result, Comp.2>-2.1)

svg('resultados.svg'), width=15, height=8)    
ggplot(result, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point(size=3) + facet_wrap(~predict, ncol=2) + ggtitle('resultado PCA+DA')
dev.off()

###other way to do it 

svg('puntitos3.svg')    
ggplot(pcaValuesSample, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point() + facet_wrap(~correcto, ncol=1) + ylim(c(-0.3,0.3))
dev.off()


svg('fig_dist.svg3')    
ggplot(foo, aes(x=Comp.1, y=Comp.2)) + geom_density2d(aes(col=site), alpha=0.3) + geom_point(aes(col=site), size=2) + facet_wrap(~site, ncol=1) + theme(legend.position='none')
dev.off()

distMetrics <- matrix(0, nrow=4, ncol=4)
rownames(distMetrics) <- c('parlamento','belen','delicias','malpica')
colnames(distMetrics) <- c('parlamento','belen','delicias','malpica')

parlamento <- subset(myData, site=='parlamento')    
belen <- subset(myData, site=='belen')    
delicias <- subset(myData, site=='delicias')    
malpica <- subset(myData, site=='malpica')    

# 1 - parlamento-belen
pairData <- rbind(parlamento,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- parlamento[sample(nrow(parlamento), sampleSize),]
sample2 <- belen[sample(nrow(belen), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['parlamento','belen'] <- conf$overall['Accuracy']
distMetrics['belen','parlamento'] <- conf$overall['Accuracy']

# 2 - malpica-parlamento
pairData <- rbind(parlamento,malpica)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- parlamento[sample(nrow(parlamento), sampleSize),]
sample2 <- malpica[sample(nrow(malpica), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['parlamento','malpica'] <- conf$overall['Accuracy']
distMetrics['malpica','parlamento'] <- conf$overall['Accuracy']

# 3 - malpica-belen
pairData <- rbind(malpica,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- malpica[sample(nrow(malpica), sampleSize),]
sample2 <- belen[sample(nrow(belen), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['malpica','belen'] <- conf$overall['Accuracy']
distMetrics['belen','malpica'] <- conf$overall['Accuracy']

# 4 - delicias-parlamento
pairData <- rbind(parlamento,delicias)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- parlamento[sample(nrow(parlamento), sampleSize),]
sample2 <- delicias[sample(nrow(delicias), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['parlamento','delicias'] <- conf$overall['Accuracy']
distMetrics['delicias','parlamento'] <- conf$overall['Accuracy']

# 5 - delicias-belen
pairData <- rbind(delicias,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- delicias[sample(nrow(delicias), sampleSize),]
sample2 <- belen[sample(nrow(belen), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['delicias','belen'] <- conf$overall['Accuracy']
distMetrics['belen','delicias'] <- conf$overall['Accuracy']

# 6 - delicias-malpica

pairData <- rbind(delicias,malpica)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- delicias[sample(nrow(delicias), sampleSize),]
sample2 <- malpica[sample(nrow(malpica), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['delicias','malpica'] <- conf$overall['Accuracy']
distMetrics['malpica','delicias'] <- conf$overall['Accuracy']
distMetrics

#######################################################CHECK WITH OTHERS POTTERY WORKSHOPS######################################################################

myData <- read.csv('dresall.csv', header=T, sep= ",")
myData= subset(myData, type %in% c("Dressel C","Dressel D","Dressel E","Dressel G","Dressel H"))
myData <- myData[,5:13]

#to select the pottery WORKSHOPS

myData=myData[myData$site=="delicias" | myData$site=="malpica" | myData$site== "belen"| myData$site== "parlamento"| myData$site== "penaflor" | myData$site== "pesebres",]  
myData$site = factor(myData$site)


##ready for PCA!!
logData <- log(myData[,1:8])
pcaResults <- princomp(logData, center=T, scale=T)
plot(pcaResults)

pcaValues <- as.data.frame(pcaResults$scores)
pcaValues$site <- myData$site
ggplot(pcaValues, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point()
ggplot(pcaValues, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point() + facet_wrap(~site,ncol=1)


sampleSize = min(count(myData,'site')$freq)

amph1 <- subset(myData, site=="delicias")
sample1 <- amph1[sample(nrow(amph1), sampleSize),]
amph2 <- subset(myData, site=="malpica")
sample2 <- amph2[sample(nrow(amph2), sampleSize),]
amph3 <- subset(myData, site=="belen")
sample3 <- amph3[sample(nrow(amph3), sampleSize),]
amph4 <- subset(myData, site=="parlamento")
sample4 <- amph4[sample(nrow(amph4), sampleSize),]
amph5 <- subset(myData, site=="pesebres")
sample5 <- amph5[sample(nrow(amph5), sampleSize),]
amph6 <- subset(myData, site=="penaflor")
sample6 <- amph6[sample(nrow(amph6), sampleSize),]

sample <- rbind(sample1,sample2)
sample <- rbind(sample, sample3)
sample <- rbind(sample, sample4)
sample <- rbind(sample, sample5)
sample <- rbind(sample, sample6)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- sample$site

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1,1,1,1,1)/6)
predQda <- predict(amphDA, pcaValuesSample)
sample$probDelicias<- predQda$posterior[,"delicias"]
sample$probDelicias<- predQda$posterior[,"malpica"]
sample$probMalpica <- predQda$posterior[,"belen"]
sample$probBelen <- predQda$posterior[,"parlamento"]
sample$probBelen <- predQda$posterior[,"pesebres"]
sample$probBelen <- predQda$posterior[,"penaflor"]
pcaValuesSample$class <- predQda$class

#CONFUSION

confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)

svg('6workshopes.svg')    
ggplot(pcaValuesSample, aes(x=Comp.1, y=Comp.2, col=site)) + geom_point() + facet_wrap(~correcto, ncol=1) + ylim(c(-0.3,0.3))
dev.off()

foo <- subset(pcaValuesSample, Comp.1>-1.2)
foo <- subset(foo, Comp.2>-2.1)

svg('fig_distseis.svg')    
ggplot(foo, aes(x=Comp.1, y=Comp.2)) + geom_density2d(aes(col=site), alpha=0.3) + geom_point(aes(col=site), size=2) + facet_wrap(~site, ncol=1) + theme(legend.position='none')
dev.off()

####TO CHECK THE DISTANCE

distMetrics <- matrix(0, nrow=6, ncol=6)
rownames(distMetrics) <- c('parlamento','belen','delicias','malpica','pesebres', 'penaflor')
colnames(distMetrics) <- c('parlamento','belen','delicias','malpica','pesebres', 'penaflor')

parlamento <- subset(myData, site=='parlamento')    
belen <- subset(myData, site=='belen')    
delicias <- subset(myData, site=='delicias')    
malpica <- subset(myData, site=='malpica')   
pesebres <- subset(myData, site=='pesebres')
penaflor <- subset(myData, site=='penaflor')

# 1 - parlamento-belen
pairData <- rbind(parlamento,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- parlamento[sample(nrow(parlamento), sampleSize),]
sample2 <- belen[sample(nrow(belen), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['parlamento','belen'] <- conf$overall['Accuracy']
distMetrics['belen','parlamento'] <- conf$overall['Accuracy']

# 2 - malpica-parlamento
pairData <- rbind(parlamento,malpica)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- parlamento[sample(nrow(parlamento), sampleSize),]
sample2 <- malpica[sample(nrow(malpica), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['parlamento','malpica'] <- conf$overall['Accuracy']
distMetrics['malpica','parlamento'] <- conf$overall['Accuracy']

# 3 - malpica-belen
pairData <- rbind(malpica,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- malpica[sample(nrow(malpica), sampleSize),]
sample2 <- belen[sample(nrow(belen), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['malpica','belen'] <- conf$overall['Accuracy']
distMetrics['belen','malpica'] <- conf$overall['Accuracy']

# 4 - delicias-parlamento
pairData <- rbind(parlamento,delicias)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- parlamento[sample(nrow(parlamento), sampleSize),]
sample2 <- delicias[sample(nrow(delicias), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['parlamento','delicias'] <- conf$overall['Accuracy']
distMetrics['delicias','parlamento'] <- conf$overall['Accuracy']

# 5 - delicias-belen
pairData <- rbind(delicias,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- delicias[sample(nrow(delicias), sampleSize),]
sample2 <- belen[sample(nrow(belen), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['delicias','belen'] <- conf$overall['Accuracy']
distMetrics['belen','delicias'] <- conf$overall['Accuracy']

# 6 - delicias-malpica

pairData <- rbind(delicias,malpica)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- delicias[sample(nrow(delicias), sampleSize),]
sample2 <- malpica[sample(nrow(malpica), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['delicias','malpica'] <- conf$overall['Accuracy']
distMetrics['malpica','delicias'] <- conf$overall['Accuracy']

# 7 pesebres- malpica

pairData <- rbind(pesebres,malpica)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- pesebres[sample(nrow(pesebres), sampleSize),]
sample2 <- malpica[sample(nrow(malpica), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['pesebres','malpica'] <- conf$overall['Accuracy']
distMetrics['malpica','pesebres'] <- conf$overall['Accuracy']

# 8 pesebres- belen

pairData <- rbind(pesebres,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- pesebres[sample(nrow(pesebres), sampleSize),]
sample2 <- belen[sample(nrow(belen), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['pesebres','belen'] <- conf$overall['Accuracy']
distMetrics['belen','pesebres'] <- conf$overall['Accuracy']

# 9 pesebres- parlamento

pairData <- rbind(pesebres,parlamento)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- pesebres[sample(nrow(pesebres), sampleSize),]
sample2 <- parlamento[sample(nrow(parlamento), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['pesebres','parlamento'] <- conf$overall['Accuracy']
distMetrics['parlamento','pesebres'] <- conf$overall['Accuracy']

# 10 pesebres- delicias

pairData <- rbind(pesebres,delicias)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- pesebres[sample(nrow(pesebres), sampleSize),]
sample2 <- delicias[sample(nrow(delicias), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['pesebres','delicias'] <- conf$overall['Accuracy']
distMetrics['delicias','pesebres'] <- conf$overall['Accuracy']


# 11 penaflor-malpica

pairData <- rbind(penaflor,malpica)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- penaflor[sample(nrow(penaflor), sampleSize),]
sample2 <- malpica[sample(nrow(malpica), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['penaflor','malpica'] <- conf$overall['Accuracy']
distMetrics['malpica','penaflor'] <- conf$overall['Accuracy']


# 12 penaflor-belen

pairData <- rbind(penaflor,belen)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- penaflor[sample(nrow(penaflor), sampleSize),]
sample2 <- belen[sample(nrow(belen), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['penaflor','belen'] <- conf$overall['Accuracy']
distMetrics['belen','penaflor'] <- conf$overall['Accuracy']

# 13 penaflor-parlamento

pairData <- rbind(penaflor,parlamento)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- penaflor[sample(nrow(penaflor), sampleSize),]
sample2 <- parlamento[sample(nrow(parlamento), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['penaflor','parlamento'] <- conf$overall['Accuracy']
distMetrics['parlamento','penaflor'] <- conf$overall['Accuracy']


# 14 penaflor- delicias

pairData <- rbind(penaflor,delicias)
sampleSize = min(count(pairData,'site')$freq)

sample1 <- penaflor[sample(nrow(penaflor), sampleSize),]
sample2 <- delicias[sample(nrow(delicias), sampleSize),]
sample <- rbind(sample1,sample2)

logDataSample <- log(sample[,1:8])
pcaResultsSample <- princomp(logDataSample, center=T, scale=T)
pcaValuesSample <- as.data.frame(pcaResultsSample$scores)
pcaValuesSample$site <- factor(sample$site)

amphDA <- qda(site~Comp.1+Comp.2, data=pcaValuesSample, prior=c(1,1)/2)
predQda <- predict(amphDA, pcaValuesSample)
pcaValuesSample$class <- predQda$class

conf <- confusionMatrix(pcaValuesSample$class, pcaValuesSample$site)
distMetrics['penaflor','delicias'] <- conf$overall['Accuracy']
distMetrics['delicias','penaflor'] <- conf$overall['Accuracy']

distMetrics











