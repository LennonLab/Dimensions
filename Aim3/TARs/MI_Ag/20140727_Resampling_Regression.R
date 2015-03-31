rm(list=ls()) #clears currernt workspace
getwd() #erases previous working directory
setwd("/Users/lennonj/Desktop/") #sets new working directory

# Generate matrix for regression
x<-matrix(rnorm(500,mean=100,sd=10))
y<-matrix(rnorm(500,mean=x,sd=10))
plot(x,y)
data<-cbind(x,y)
colnames(data)<-c("x","y")
data<-as.data.frame(data)

# Jack-knifing
R<-999 # number of interations
out<-matrix(data=NA,nrow=R,ncol=2) # create output matrix for slope and intercept

# Call regression function to get intercepts and slopes, respectively
reg.parms<-function(x,y) {
	reg<-lm(y~x,data=jack)
	matrix.out<-cbind(summary(reg)$coefficients[1],summary(reg)$coefficients[2])
	return(matrix.out)
}	
	
# For-loop to populate out matrix
# samples n obs (i.e., "size") from original data with replacement
for (i in 1:R) {
	jack<-data[sample(1:nrow(data),size=400,replace=T),] 
	out[i,1:2] <- reg.parms(jack$x, jack$y)
} 

# Distributions from jackknifed data
colnames(out)<-c("intercepts","slopes")
boxplot(out[,1],ylab="intercept")
boxplot(out[,2],ylab="slope")

# Confidence intervals 
intercept_CI<-quantile(out[,1],c(0.025,0.975)) 
slope_CI<-quantile(out[,2],c(0.025,0.975)) 

