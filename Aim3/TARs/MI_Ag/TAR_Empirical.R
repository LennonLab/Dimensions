## Set WD; /Tar Project File/Empirical Data
rm(list=ls()) #clears currernt workspace
getwd() #erases previous working directory
setwd("/Users/lennonj/GitHub/TARs/Empirical_TAR")
#setwd("/Users/Kayla I Miller/My Box Files/TAR Project File/Empirical Data") # Change to GitHub directory
# test

## Loading Packages 
library(vegan)
library(gmt)
library(som)
library(boot)

## Bringing in matrix from mothur ##
source("Functions/DiversityFunctions.r") 
shared = "Input/AG_TAR.july.shared"
design = "Input/AGTAR.design"

# Import Site by OTU Matrix
# Option: Change OTU Cut-off (0.03, 0.05, or 0.01)
AGTAR03.raw <- t(read.otu(shared, "0.03")) 

# Bring in design file
design <- read.delim(design, header=T, row.names=1)

# removesOTUs found once (or less)
AGTAR03 <- AGTAR03.raw[rowSums(AGTAR03.raw) > 1,]

#AGTAR01 (99%) - 134,590 -> 72,961
#AGTAR03 (97%) - 33,055 -> 24,716
#AGTAR05 (95%) - 17,433 -> 13,960 

####
# Calculate Presence Absence
AGTAR03.PA <- (AGTAR03 > 0)*1

####
# Calculate Relative Abundance
AGTAR03.REL <- AGTAR03
for(i in 1:ncol(AGTAR03)){
    AGTAR03.REL[,i] = AGTAR03[,i]/sum(AGTAR03[,i])
    } 
	
# Taking log of rel.abundance
# getting error, but maybe still works? should confirm
logAGTAR03.REL <- decostand(AGTAR03.REL, method="log",logbase=10) 

##################################################################
# Making a Dormant datafile 

## Changing the matrix, so make a copy: 
DormMatPA <-  AGTAR03.PA

# Counter for going through columns (i)
# These are the columns in my dataset that correspond to DNA samples
DNA <- c(1:11, 23:33, 45:55)

# Getting the number of rows for a counter & setting a counter to 0 

j <- 0 

for (i in DNA){ ## i is the col number for DNA
	j = i + 11 ## the corresponding cDNA column is 11 cols away
	for (y in 1:nrow(DormMatPA)) {  ## for each row in the matrix
		if (DormMatPA[y,j] == 1) {  ## if the row in the RNA col = 1 (active)
			DormMatPA[y,i] <- (DormMatPA[y,i]) * 0  ## then multiply the DNA column by 0 to make remove active OTU from dormant
}}}

# compare head(AGTAR03.PA) with head(DormMatPA); DNA in latter should == 0

### Making a Dormant datafile 
## !! Relative Abundance !! ##

DormMatREL <-  AGTAR03.REL
j <- 0

for (i in DNA){ ## i is the col number for DNA
	j = i + 11 ## the corresponding cDNA column is 11 cols away
	for (y in 1:nrow(DormMatREL)) {  ## for each row in the matrix
		if (DormMatREL[y,j] > 0) {  ## if the row in the RNA col = 1 (active)
			DormMatREL[y,i] <- (DormMatREL[y,i]) * 0  ## then multiply the DNA column by 0 to make remove active OTU from dormant
}}}

# compare head(AGTAR03.REL) with head(DormMatREL); DNA in latter should == 0

#######################################################################
## Geographic Distance 
## Adapted From Zinger
## Distances are in km 

###### Seems to be a potential error in geogrpahic distances
###### I thought smaller distance was 1-10 cm
###### min(dist.ag.vec) = 0.001 == 100 cm
###### noticed this when trying to fix x-axis labels in DD plots

## AG COORDS
coords.ag <- read.table("Input/coords_half.txt", header=TRUE) 
# possibly add more decimals to lat/longs to avoid double zeroes?
rownames = list(coords.ag$cid)

dists.ag = matrix(0, nrow(coords.ag), nrow(coords.ag), dimnames=list(rownames(rownames),rownames(rownames)))

for(j in 1:(nrow(coords.ag)-1)) for(i in (j+1):nrow(coords.ag)){ #for each sample pairwise comparison:
  dists.ag[i,j] = (geodist(as.numeric(as.vector(coords.ag[i,"LAT"])), as.numeric(as.vector(coords.ag[i,"LONG"])), as.numeric(as.vector(coords.ag[j,"LAT"])), as.numeric(as.vector(coords.ag[j,"LONG"])), units="km"))
  }

####
## Extracting values from matricies into vector

dist.ag.vec <- NULL
j <- 1 # j = col
for (i in 2:33)   {  # i = row
	dist.ag.vec <- append(dist.ag.vec, dists.ag[i:33, j]) 
	j = j + 1
	i = i+1
}

dist.ag.vec <- dist.ag.vec + 0.0000001 

lndist.ag <- log(dist.ag.vec) # log transforming 
 
########################################################################
## Community Distance ##
# Options: Different distance methods; using PA vs. RelAbund 
# Vegdist returns dissimiliarities 
# Options for next line: Binary "TRUE" = Sorenson; "FALSE" = Bray
# Vegdist wants rows to be OTUs, so need to transpose

pre_dormmat.bray = as.matrix(vegdist(t(DormMatREL),method="bray",binary = FALSE))

#dormmat.bray <-log(1-(pre_dormmat.bray)) # use to ln-transform similarity
dormmat.bray <-(1-(pre_dormmat.bray)) # 1 - dissimilarity ==> give similarity
dormmat.bray.DNA <- dormmat.bray[c(1:11, 23:33, 45:55),c(1:11, 23:33, 45:55)]
dormmat.bray.RNA <- dormmat.bray[c(12:22, 34:44, 56:66),c(12:22, 34:44, 56:66)]

## total 
pre_total.bray = as.matrix(vegdist(t(AGTAR03.REL),method="bray",binary = FALSE))

total.bray <-(1-(pre_total.bray)) # 1 - dissimilarity
total.bray <- total.bray[c(1:11, 23:33, 45:55),c(1:11, 23:33, 45:55)]

########################################################################
## Estimating Richness of Dormant and Active Communities
## This is JTL approach. Need to think about it some more
## Very low number of seqs in some dorm samples after "subtracting" active

Dorm.Seqs <-  AGTAR03 # Make raw seq matrix for esimating richness
j <- 0

for (i in DNA){ ## i is the col number for DNA
	j = i + 11 ## the corresponding cDNA column is 11 cols away
	for (y in 1:nrow(Dorm.Seqs)) {  ## for each row in the matrix
		if (Dorm.Seqs[y,j] > 0) {  ## if the row in the RNA col = 1 (active)
			Dorm.Seqs[y,i] <- (Dorm.Seqs[y,i]) * 0  ## then multiply the DNA column by 0 to make remove active OTU from dormant
}}}

#### Percent sequences in dormant and active pools
sum.total.seqs<-as.matrix(t(colSums(AGTAR03))) # total # of seqs after removing global singletons
DNA.seqs<-AGTAR03[,c(1:11, 23:33, 45:55)] # pull DNA seqs out
RNA.seqs<-AGTAR03[,c(12:22, 34:44, 56:66)] # pull RNA seqs out
Dorm.seqs.PA<-(Dorm.Seqs > 0)*1 # presence-absence matrix of Dorm matrix
Dorm.MAT<-Dorm.seqs.PA[,c(1:11, 23:33, 45:55)] # pull out DNA PA 
Act.MAT<-Dorm.seqs.PA[,c(12:22, 34:44, 56:66)] # pull out RNA PA 
Dorm.abund<-DNA.seqs * Dorm.MAT # multiply dorm matrix by seq matrix
Act.abund<-DNA.seqs * Act.MAT # multiply active matrix by seq matrix
sum.DNA.seqs<-sum.total.seqs[c(1:11, 23:33, 45:55)]
sum.RNA.seqs<-sum.total.seqs[c(12:22, 34:44, 56:66)]
sum.Dorm.seqs<-colSums(Dorm.abund)
sum.Act.seqs<-colSums(Act.abund)
per.Dorm.seqs<-(sum.Dorm.seqs/sum.DNA.seqs)*100
per.Act.seqs<-(sum.Act.seqs/sum.DNA.seqs)*100
per.seqs.means<-c(mean(per.Dorm.seqs),(mean(per.Act.seqs)))
sem <- function(x, ...){sd(x)/sqrt(length(x))}
per.seqs.sem<-c(sem(per.Dorm.seqs),(sem(per.Act.seqs)))
per.seqs.x<-c("dormant","active")
per.seqs.t<-t.test(per.Dorm.seqs,per.Act.seqs)

# Barplot for percent dormant sequences
png(filename="./Plots/Dorm_Seqs.png",width=1200, height=1200, res=96*2) # this says following plot will be saved to folder
par(bg="white",mar=c(7,7,1,1),lwd=3)
per.seqs.plot <- barplot(per.seqs.means, names.arg=per.seqs.x, ylim=c(0,100), 
col=c("yellow","slategray1"),yaxt="n", cex.names=2)
arrows(x0 = per.seqs.plot, y0 = per.seqs.means, y1 = per.seqs.means - per.seqs.sem, angle = 90, 
length=0.25, lwd = 2) # lower error bars
arrows(x0 = per.seqs.plot, y0 = per.seqs.means, y1 = per.seqs.means + per.seqs.sem, angle = 90, 
length=0.25, lwd = 2) # upper error bars
axis(side = 2, labels=T, lwd.ticks=3, las=2, lwd=3, cex=2, cex.axis=2)
mtext("Percent of Sequences (%)", side = 2, cex = 2, line = 5.5)
dev.off() # this writes plot to folder
graphics.off() # may need to use this so you can plot to screen

#### OTUs in dormant and active pools
OTU.dorm<-colSums(Dorm.MAT)
OTU.act<-colSums(Act.MAT)
OTU.tot<-OTU.dorm+OTU.act
OTU.per.dorm<-(OTU.dorm/OTU.tot)*100
OTU.means<-c(mean(OTU.dorm),mean(OTU.act))
OTU.sem<-c(sem(OTU.dorm),sem(OTU.act))
OTU.x<-c("dormant","active")
OTU.t<-t.test(OTU.dorm,OTU.act)

# Barplot for OTU in dormant and active
png(filename="./Plots/Dorm_OTUs.png",width=1200, height=1200, res=96*2) # this says following plot will be saved to folder
par(bg="white",mar=c(7,7,1,1),lwd=3)
OTU.plot <- barplot(OTU.means, names.arg=OTU.x, ylim=c(0,4000), 
col=c("yellow","slategray1"),yaxt="n", cex.names=2)
arrows(x0 = OTU.plot, y0 = OTU.means, y1 = OTU.means - OTU.sem, angle = 90, 
length=0.25, lwd = 3) # lower error bars
arrows(x0 = OTU.plot, y0 = OTU.means, y1 = OTU.means + OTU.sem, angle = 90, 
length=0.25, lwd = 3) # upper error bars
axis(side = 2, labels=T, lwd.ticks=3, las=2, lwd=3, cex=2, cex.axis=2)
mtext("Taxa Richness", side = 2, cex = 2, line = 5.5)
dev.off() # this writes plot to folder
graphics.off() # may need to use this so you can plot to screen

#### Pielou's Evenness in dormant and active pools
even.dorm<-(diversity(t(Dorm.abund)))/log(OTU.dorm)
even.act<-(diversity(t(Act.abund)))/log(OTU.act)
even.means<-c(mean(even.dorm),mean(even.act))
even.sem<-c(sem(even.dorm),sem(even.act))
even.x<-c("dormant","active")
even.t<-t.test(even.dorm,even.act)

#### Evar Evenness in dormant and active pools
even.dorm<-Evar(t(Dorm.abund))
even.act<-Evar(t(Act.abund))
even.means<-c(mean(even.dorm),mean(even.act))
even.sem<-c(sem(even.dorm),sem(even.act))
even.x<-c("dormant","active")
even.t<-t.test(even.dorm,even.act)

# Barplot for evenness in dormant and active
png(filename="./Plots/Dorm_even.png",width=1200, height=1200, res=96*2) # this says following plot will be saved to folder
par(bg="white",mar=c(7,7,1,1),lwd=3)
even.plot <- barplot(even.means, names.arg=even.x, ylim=c(0,1), 
col=c("yellow","slategray1"),yaxt="n", cex.names=2)
arrows(x0 = even.plot, y0 = even.means, y1 = even.means - even.sem, angle = 90, 
length=0.25, lwd = 3) # lower error bars
arrows(x0 = even.plot, y0 = even.means, y1 = even.means + even.sem, angle = 90, 
length=0.25, lwd = 3) # upper error bars
axis(side = 2, labels=T, lwd.ticks=3, las=2, lwd=3, cex=2, cex.axis=2)
mtext("Taxa Evenness (Evar)", side = 2, cex = 2, line = 5.5)
dev.off() # this writes plot to folder
graphics.off() # may need to use this so you can plot to screen

########################################################################
## Extracting values from matricies into vectors
## And calculating SIMILARITIES 

dormmat.bray.DNA.dis <- NULL
j <- 1 # j = col
for (i in 2:33)   {  # i = row; counter for rows; starting with 2 since there are 0s along diagonal
	dormmat.bray.DNA.dis<- append(dormmat.bray.DNA.dis, dormmat.bray.DNA[i:33, j]) 
	j = j + 1
	i = i+1
}

## same as before for RNA 
dormmat.bray.RNA.dis <- NULL
j <- 1 # j = col
for (i in 2:33)   {  # i = row; counter for rows; starting with 2 since there are 0s along diagonal
	dormmat.bray.RNA.dis<- append(dormmat.bray.RNA.dis, dormmat.bray.RNA[i:33, j]) 
	j = j + 1
	i = i+1
}

## total
total.bray.dis <- NULL
j <- 1 # j = col
for (i in 2:33)   {  # i = row; counter for rows; starting with 2 since there are 0s along diagonal
	total.bray.dis<- append(total.bray.dis, total.bray[i:33, j]) 
	j = j + 1
	i = i+1
}

#######################################################3################
## Calculating linear regressions

dormmat.bray.RNA.line <- lm(dormmat.bray.RNA.dis~lndist.ag)
dormmat.bray.DNA.line <- lm(dormmat.bray.DNA.dis~lndist.ag)
total.bray.line <- lm(total.bray.dis~lndist.ag)

dormmat.bray.DNA.dis <- dormmat.bray.DNA.dis
dormmat.bray.RNA.dis <- dormmat.bray.RNA.dis
total.bray.dis <- total.bray.dis

########################################################################
## Plotting

png(filename="./Plots/OTU_Decay.png",width=1200, height=1200, res=96*2) # this says following plot will be saved to folder

par(bg="white",mar=c(5,5,1,1)) # make plot background white # mar: bottom,left,top,right
plot(dist.ag.vec, dormmat.bray.DNA.dis, log="x", xaxt="n",yaxt="n",xlab="Distance (km)", ylab= "OTU Similarity (Bray-Curtis)", cex.lab=1.5, pch=22, col="black", bg="yellow", lwd=1.8, cex=2.4, xlim=c(0.00005, 200),ylim=c(0,1),las=2, bty="n",box(lwd=2)) # Dormant points
axis(side=2,lwd=2,las=2,cex.axis=1.5)
#ticks <- seq(-5, 2, by=1)
ticks <- seq(-4, 2, by=2) # creates tick range and increment
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i)))) # generates label exponents
#axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100), lwd=2,cex.axis=1.5,labels=labels)
axis(1, at=c(0.0001, 0.01, 1, 100), lwd=2,cex.axis=1.5,labels=labels) # identifies labe values

# Following plots regression lines for a range of X
dormmat.line <- lm(dormmat.bray.DNA.dis~dist.ag.vec) # function to generate predicteds 
dorm.reg.range<-c(0.0001,150) # range to plot predicted
dorm.sim<-predict(dormmat.line,newdata=data.frame(dist.ag.vec=dorm.reg.range)) # predicted to plot
lines(x=dorm.reg.range,y=dorm.sim,lty=2,lwd=3) # plot

points(dist.ag.vec, dormmat.bray.RNA.dis, pch=22, col="black",bg="slategray1",lwd=1.8, cex=2.4) # Active points
active.line <- lm(dormmat.bray.RNA.dis~dist.ag.vec) # function to generate predicteds 
active.reg.range<-c(0.0001,150) # range to plot predicted
active.sim<-predict(active.line,newdata=data.frame(dist.ag.vec=active.reg.range)) # predicted to plot
lines(x=active.reg.range,y=active.sim,lty=2,lwd=3) # plot

#legend("topright", inset=0.05,c("active","dormant"),lty=1, lwd= 5, col=c("yellow","slategray1"),bty="n")

legend("topright", inset=0.05,c("active","dormant"),pch=22, cex=1.5, col=c("black","black"), pt.bg=c("slategray1","yellow"),pt.cex=2.5, lwd= 1.8, ,bty="n",lty=NA)

dev.off() # this writes plot to folder
graphics.off() # may need to use this so you can plot to screen

#points(dist.ag.vec, total.bray.dis, pch=22, col="black",bg="purple",lwd=1.8,cex=2.4) # Total points

########################################################################
## Bootstrapping

R<-999 # number of interations
boot.out<-matrix(data=NA,nrow=R,ncol=4) # create output matrix for slope and intercept

dorm.boot<-data.frame(lndist.ag,dormmat.bray.DNA.dis)
colnames(dorm.boot)<-c("x","y")
active.boot<-data.frame(lndist.ag,dormmat.bray.RNA.dis)
colnames(active.boot)<-c("x","y")

# Call regression function to get intercepts and slopes, respectively
reg.parms<-function(x,y) {
	reg<-lm(y~x,data=jack)
	matrix.out<-cbind(summary(reg)$coefficients[1],summary(reg)$coefficients[2])
	return(matrix.out)
}	

# For-loop to populate out matrix
# samples n obs (i.e., "size") from original data with replacement

for (i in 1:R) {
	jack<-dorm.boot[sample(1:nrow(dorm.boot),size=500,replace=T),] 
	boot.out[i,1:2] <- reg.parms(jack$x, jack$y)
} 

for (i in 1:R) {
	jack<-active.boot[sample(1:nrow(active.boot),size=500,replace=T),] 
	boot.out[i,3:4] <- reg.parms(jack$x, jack$y)
} 

colnames(boot.out)<-c("dormant","dormant","active","active")

ints <-boot.out[,c(1,3)]
slopes <- boot.out[,c(2,4)]

## Plotting Box and Wisker plots 
par(mfrow=c(1,2),oma=c(0,0,0,0))
boxplot(ints, main = "intercepts", ylab = "Similarity",las=2)
boxplot(slopes, main = "slopes", ylab = "Decay (Similarity / ln Km)",las=2)

## Standard error calculations
se <- function(x) sqrt(var(x)/length(x))
se(ints)
se(slopes)

## 95% confidence intervals 
dorm.slope.LCI<-quantile(slopes[,1],0.025)
dorm.slope.UCI<-quantile(slopes[,1],0.975)
active.slope.LCI<-quantile(slopes[,2],0.025)
active.slope.UCI<-quantile(slopes[,2],0.975)

# Significant test: is the max of LCI less than or equal to min of UCI?
# If "FALSE", then confidence intervals do not overlap
max(dorm.slope.LCI,active.slope.LCI) <= min(dorm.slope.UCI,active.slope.UCI)

dorm.int.LCI<-quantile(ints[,1],0.025)
dorm.int.UCI<-quantile(ints[,1],0.975)
active.int.LCI<-quantile(ints[,2],0.025)
active.int.UCI<-quantile(ints[,2],0.975)

max(dorm.int.LCI,active.int.LCI) <= min(dorm.int.UCI,active.int.UCI)

### ---> JTL hasn't done much with code below here
### ---> Seems like there's some duplicaiton of big chunks
### ---> Needs some curating

########################################################################
# Principal Coordinates Analysis
# Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
# eig=TRUE returns eigenvalues; k = # of dimensions to calculate

##Input is a distance matrix - options here, too
AGTAR_pcoa <- cmdscale(pre_dormmat.bray, k=3, eig=TRUE, add=FALSE)

# Percent Variance Explained Using PCoA (Axis 1,2,3)
explainvar1 <- round(AGTAR_pcoa$eig[1]/sum(AGTAR_pcoa$eig)*100,2) 
explainvar2 <- round(AGTAR_pcoa$eig[2]/sum(AGTAR_pcoa$eig)*100,2)
explainvar3 <- round(AGTAR_pcoa$eig[3]/sum(AGTAR_pcoa$eig)*100,2)

pcoap <- merge(as.data.frame(AGTAR_pcoa$points),design,by=0,all.x=T)[,-1]
rownames(pcoap) <- rownames(AGTAR_pcoa$points)

# Plot Parameters
par(mfrow=c(1,1), mar=c(5,5,1,1)) 
layout(rbind(1, 2), height=c(7, 1)) 
y.dim <- c(min(pcoap$V2)+min(pcoap$V2)*0.2,max(pcoap$V2)+max(pcoap$V2)*0.2)
x.dim <- c(min(pcoap$V1)+min(pcoap$V1)*0.2,max(pcoap$V1)+max(pcoap$V1)*0.2)
  
# Initiate Plot
plot(pcoap$V1, pcoap$V2, xlab=paste("PCoA Axis 1 (",explainvar1, "%)", sep="")
    , ylab=paste("PCoA Axis 2 (",explainvar2, "%)", sep=""), 
    xlim=x.dim,ylim= y.dim, pch=16, cex=2.0, type="n",xaxt="n",
    yaxt="n", cex.lab=1.5, cex.axis=1.2, )  
axis(side=1, las=1)   
axis(side=2, las=1)    
abline(h=0, lty="dotted")  
abline(v=0, lty="dotted")

## Site variables
site.shape <- rep(NA, dim(pcoap)[1])
    for (i in 1:length(site.shape)){
      if (pcoap$site[i] == "EL"){site.shape[i] = 24}
	else if (pcoap$site[i] == "KBS"){site.shape[i] = 21}
	else if (pcoap$site[i] == "RF"){site.shape[i] = 22}
      else {site.shape[i] = 3}
      }

#nucleic acid variables 
natype.color <- rep(NA, dim(pcoap)[1])
    for (i in 1:length(natype.color)){
      if (pcoap$natype[i] == "DNA") {natype.color[i] = "blue"}
      else if (pcoap$natype[i] == "RNA") {natype.color[i] = "red"}
      else {natype.color[i] = "white"}
      } 

points(pcoap$V1, pcoap$V2, pch=site.shape, cex=2.0, bg=natype.color, lwd=1)   
 
# Legend 
box(lwd=2)
par(mar=c(0, 3, 0, 0))
plot.new()
legend("center", c(paste("East Lansing", sep=""), 
 paste("KBS"), 
 paste("Russ Forest"),
 paste ("Dormant"),
 paste ("Active")), 
pt.lwd=2, col="black", pt.bg=c("black", "black", "black", "blue", "red"), pch=c(24,21,22,23,23), bty='n', ncol=2, cex=1.5, pt.cex=2)  

########################################################################
## PERMANOVA
## Input is a distance matrix, same options here
## Method = tell it which method you used in vegdist

Adonis_Dorm <- adonis(dormmat.bray ~ design$site*design$natype, method="bray", 
    permutations=1000)
Adonis_Dorm

########################################################################
## Calculating % dormancy
## Adapted Code from Mario Muscarella 

## For counter (i); columns that contain DNA samples 
DNA <- c(1:11, 23:33, 45:55)

per.d  <- rep(NA, 33)
total  <- rep(NA, 33)
active <- rep(NA, 33)
site <- rep(NA,33)

j <- 0
y <- 0

for (i in DNA){

		j = i + 11  ## Corresponding RNA sample is 11 away from DNA 
		y = y+1

		X.tmp1 <- AGTAR03.PA[, c(i,j)] # Making a temporary matrix with a paired sample (DNA & RNA)
		X.tmp2 <- X.tmp1[which(rowSums(X.tmp1) != 0),]  ## getting rid of 0 rows in either (should be both?) DNA/RNA
		tot <- dim(X.tmp2)[1] ## counting total number of rows 
		total[y] <- tot
		act <- sum(X.tmp2[,2]) ## counting number of rows in column 2 (which should be RNA, based on alphabet)
		active[y] <- act
		per.d[y] <- 1 - (act/tot)
		site[y] <- row.names(t(AGTAR03.PA[, c(i,j)])) ## Adding one of the site labels as it iterates
}
# Because my mind prefers to look at it like this
perc <- per.d * 100 
# adding all to one data frame for ease 
PerDorm <- data.frame(site, total, active, per.d, perc)

# some stats 
mean(perc)
median(perc)
range(perc)

########################################################################
##### Environmental Heterogeneity ### 
### Euclidean Distance

## Bringing in data
env_ag_full <- read.table("Input/env_ag.txt", header=TRUE)
env_for_full <- read.table("Input/env_for.txt", header=TRUE)

## removing columns (location ID)
env_ag <- env_ag_full[,-1]
env_for <- env_for_full[,-1]

# rownames for later
rownames2 = list(env_ag_full$lid)
rownames.for = list(env_for_full$lid)

## Standardizing data
## SOM package; mean = 0; variance = 1
## Gives negative values, errors in BC distances

env_ag_norm <- normalize(env_ag, byrow=FALSE)
env_for_norm <- normalize(env_for, byrow=FALSE)

## Another option 
## From Vegan - decostand; range = numbers between 0-1 
## Can use BC/Jaccard then 

env_for_01 <- decostand(env_for, method="range")
env_ag_01 <- decostand(env_ag, method="range")

## Bringing in Geographic Distances 
## Forest Coords

coords.for <- read.table("Input/coords_for.txt", header=TRUE) 
rownames.for = list(coords.for$cid)

dists.for = matrix(0, nrow(coords.for), nrow(coords.for), dimnames=list(rownames(rownames.for),rownames(rownames.for)))

for(j in 1:(nrow(coords.for)-1)) for(i in (j+1):nrow(coords.for)){ #for each sample pairwise comparison:
  dists.for[i,j] = (geodist(as.numeric(as.vector(coords.for[i,"LAT"])), as.numeric(as.vector(coords.for[i,"LONG"])), as.numeric(as.vector(coords.for[j,"LAT"])), as.numeric(as.vector(coords.for[j,"LONG"])), units="km"))
  }

####
## Extracting values from matricies 

dist.for.vec <- NULL
j <- 1 # j = col
for (i in 2:33)   {  # i = row
	dist.for.vec <- append(dist.for.vec, dists.for[i:33, j]) 
	j = j + 1
	i = i+1
}
##add 0.001 so no values are zero
dist.for.vec <- dist.for.vec + 0.001

lndist.for <- log(dist.for.vec)

## Dissimilarity matrix 
## AG
env.ag.bray.matrix = matrix(0, nrow(env_ag), nrow(env_ag), dimnames=list(rownames(rownames2),rownames(rownames2)))

env.ag.bray.matrix = as.matrix((vegdist(env_ag_01, method="bray")))

env.ag.bray.vec <- NULL
j <- 1 # j = col
for (i in 2:33)   {  # i = row
	env.ag.bray.vec <- append(env.ag.bray.vec, env.ag.bray.matrix[i:33, j]) 
	j = j + 1
	i = i+1
}

env.ag.bray.dis <- 1-(env.ag.bray.vec) 

lmenv.bray.ag <- lm(env.ag.bray.dis~lndist.ag)

## FOR
env.for.bray.matrix = matrix(0, nrow(env_for), nrow(env_for), dimnames=list(rownames(rownames.for),rownames(rownames.for)))

env.for.bray.matrix = as.matrix((vegdist(env_for_01, method="bray")))


env.for.bray.vec <- NULL
j <- 1 # j = col
for (i in 2:33)   {  # i = row
	env.for.bray.vec <- append(env.for.bray.vec, env.for.bray.matrix[i:33, j]) 
	j = j + 1
	i = i+1
}
env.for.bray.dis <- 1-(env.for.bray.vec) 

lmenv.bray.for <- lm(env.for.bray.dis~lndist.for)

## Plotting 

par(mfrow=c(1,2))
plot(lndist.ag, env.ag.bray.dis, main="Row Crop Agriculture", xlab="Geographic Distance (ln(km))", ylab="Bray Curtis Dissimilarity  ", pch=16, col="brown")
abline(lmenv.bray.ag, col="brown")

plot(lndist.for, env.for.bray.dis, main="Old Growth Forests", xlab="Geographic Distance (ln(km))", ylab="Bray Curtis Dissimilarity  ", pch=18, col="green")
abline(lmenv.bray.for, col="green")
title(main="Distance Decay of Environmental Heterogeneity")
legend("topright", inset=0.05,c("active","dormant"),lty=1,col=c("red","blue"),bty="n")


### Jay's plotting for ESA

png(filename="./Plots/Env_Decay.png",width=1200, height=1200, res=96*2) # this says following plot will be saved to folder
par(bg="white",mar=c(5,5,1,1)) # make plot background white # mar: bottom,left,top,right
plot(dist.ag.vec, env.ag.bray.dis, log="x", xaxt="n",yaxt="n",xlab="Distance (km)", ylab= "Environmental Similarity (Bray-Curtis)", cex.lab=1.5, pch=22, col="black", bg="white", lwd=1.8, cex=2.4, xlim=c(0.00005, 200),ylim=c(0,1),las=2, bty="n",box(lwd=2)) # Dormant points
axis(side=2,lwd=2,las=2,cex.axis=1.5)
#ticks <- seq(-5, 2, by=1)
ticks <- seq(-4, 2, by=2) # creates tick range and increment
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i)))) # generates label exponents
#axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100), lwd=2,cex.axis=1.5,labels=labels)
axis(1, at=c(0.0001, 0.01, 1, 100), lwd=2,cex.axis=1.5,labels=labels) # identifies labe values

# Following plots regression lines for a range of X
env.line <- lm(env.ag.bray.dis~dist.ag.vec) # function to generate predicteds
env.reg.range<-c(0.0001,150) # range to plot predicted
env.sim<-predict(env.line,newdata=data.frame(dist.ag.vec=env.reg.range)) # predicted to plot
lines(x=env.reg.range,y=env.sim,lty=2,lwd=3) # plot

dev.off() # this writes plot to folder
graphics.off() # may need to use this so you can plot to screen

## Mantel tests
mant_dorm<- pre_dormmat.bray[c(1:11, 23:33, 45:55),c(1:11, 23:33, 45:55)]
mant_act <- pre_dormmat.bray[c(12:22, 34:44, 56:66),c(12:22, 34:44, 56:66)]
mant_env<- env.ag.bray.matrix
mant_geo<-log(dists.ag + 0.0000001) 

## corr between geo distance and OTU similarity
mantel(mant_dorm,mant_env) # dorm vs. env: 0.4294, 0.001
mantel(mant_act,mant_env) # act vs. env: 0.6609, 0.001

## corr between geo distance and OTU similarity after removing environment effect
mantel.partial(mant_dorm,mant_geo, mant_env, method = "pearson", permutations = 999) # r: 0.56, 0.001 
mantel.partial(mant_act,mant_geo, mant_env, method = "pearson", permutations = 999) # r: 0.44, 0.001 







####### Kayla's Diversity functions

## Simple OTU counting
## Starting a new data.frame for output 

sample <- row.names(t(AGTAR03.PA))
riches <- data.frame(sample) 

REL.dorm.riches.99 <- rep(NA,66) ## start a vector for each OTU matrix
y <- 0 # Start counter for output

for (i in 1:66){
		y <- y+1
		tmp1 <- Dorm01.REL[,i] # Going through each colum of the OTU matrix
		tmp2 <- which(tmp1 != 0)  ## getting rid of OTUs with 0s
		rich <- length(tmp2) ## counting total number of rows/OTUs present
		REL.dorm.riches.99[y] <- rich ## Adding that value to output vector
}

riches = cbind(riches, REL.dorm.riches.99) ## Adding each vector to data frame; each vector = 1 OTU matrix

write.csv(riches, "Output/riches.csv") ##outputting data

## Vegdist calculations - simpson/shannon

td01.REL <- t(Dorm01.REL)
dorm.01.REL <- diversity(td01.REL, index="simpson") 
simpson <- cbind(simpson, dorm.01.REL)

td01.PA <- t(Dorm01.PA)
dorm.01.PA <- diversity(td01.PA, index="simpson") 
simpson <- cbind(simpson, dorm.01.PA)

write.csv(shannon, "Output/shannon.csv")


## Plotting 

riches <- read.csv("Output/rich_sum.csv", header=TRUE)


## Subsetting data

tot <- riches[which(riches$act == "total"),] 
act <- riches[which(riches$act == "active"),] 
dorm <- riches[which(riches$act == "dormant"),] 


par(mfrow=c(1,3))


boxplot(act$X95, dorm$X95, tot$X95, xaxt="n", xlab="Activity", ylab=" # of OTUS", main="95% OTU Cutoff")

axis(
  1, # puts the axis at the bottom
  at=1:3, # labels will be placed in the 3 categories
  labels=c("Active", "Dormant", "Total"), # labels 
  lwd=0, # width of the long axis line is zero, makes invisible
  lwd.ticks=0, # width of the etick lines also zero, makes them invisible
  cex.axis=1, # offset from the axis of the labels
  mgp=c(0,0,0), # middle zero controls distance of labels from axis
  )

boxplot(act$X97, dorm$X97, tot$X97, xaxt="n", xlab="Activity", ylab=" # of OTUS", main="97% OTU Cutoff")

axis(
  1, # puts the axis at the bottom
  at=1:3, # labels will be placed in the 3 categories
  labels=c("Active", "Dormant", "Total"), # labels 
  lwd=0, # width of the long axis line is zero, makes invisible
  lwd.ticks=0, # width of the etick lines also zero, makes them invisible
  cex.axis=1, # offset from the axis of the labels
  mgp=c(0,0,0), # middle zero controls distance of labels from axis
  )


boxplot(act$X99, dorm$X99, tot$X99, xaxt="n", xlab="Activity", ylab=" # of OTUS", main="99% OTU Cutoff")

axis(
  1, # puts the axis at the bottom
  at=1:3, # labels will be placed in the 3 categories
  labels=c("Active", "Dormant", "Total"), # labels 
  lwd=0, # width of the long axis line is zero, makes invisible
  lwd.ticks=0, # width of the etick lines also zero, makes them invisible
  cex.axis=1, # offset from the axis of the labels
  mgp=c(0,0,0), # middle zero controls distance of labels from axis
  )