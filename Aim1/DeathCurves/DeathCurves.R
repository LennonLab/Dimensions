## Retrieve and Set Your Working Directory
rm(list = ls())
getwd()
setwd("~/GitHub/Dimensions/Aim1/DeathCurves/")
#setwd("/Files/NDongoing/DimensionsBiodiversity/persistence_curves")

## Load Data
obs <- read.csv("longtermdormancy_20150425.csv", header = TRUE, stringsAsFactors = FALSE)

## Cleaning Up Data
#obs[1876:1879,1]="14-Feb-13"
#obs[1876:1879,2]="19-Dec-13"
#obs[1432:1439,2]="25-Sep-13"
#obs[1495:1498,2]="4-Oct-13"
#obs=obs[obs$X!="too many to count",]
#obs=obs[obs$Strain!="",]
obs <- obs[obs$Dilution!="50X",]
#obs=obs[obs$Colonies!="fungus",]
#obs=obs[obs$Colonies!="dropped due to fungal contamination",]
#obs=obs[,-8]
#obs$Dilution[obs$Dilution=="-1"]=1
#obs$Dilution[obs$Dilution=="-1.7"]=1.7

## Estimating CFUs
## Adding 1 to deal with log(0) observations --> should we just remove instead?

obs$Abund <- as.numeric(obs$Colonies) * 10 ^ as.numeric(obs$Dilution) + 1

# likelihood function for linear fit to log transformed data
fitLogLinearDecay<-function(p,N,time){
	k=p[1]
	N0=p[2]
	sd=exp(p[3])
	
	
	Nt=N0-k*time
	
	-sum(dnorm(N,Nt,sd,log=TRUE))
}

# likelihood function for quadratic fit to log transformed data
fitLogQuadDecay<-function(p,N,time){
	k2=p[1]
	k=p[2]
	N0=p[3]
	sd=exp(p[4])

	Nt=N0-k2*time^2-k*time
	
	-sum(dnorm(N,Nt,sd,log=TRUE))
}


strains=sort(unique(obs$Strain))
#only fit strains with more than 10 observations
strains=strains[table(obs$Strain)>10]
obs=obs[obs$Strain%in%strains,]

#pdf('decayFits.pdf')

# matrix for storing model fit output
summ=matrix(NA,length(strains)*max(obs$Rep),10)
counter=1

for(i in 1:length(strains)){
	strainObs=obs[obs$Strain==strains[i],]
	
	reps=unique(strainObs$Rep)
	for(j in 1:length(reps)){
		
		repObs=strainObs[strainObs$Rep==reps[j],]

		if(nrow(repObs)>3){
			
		start=repObs[1,1]
		
		time=(as.numeric(strptime(repObs$Firstread_date,format="%d-%b-%y",tz="EST"))-as.numeric(strptime(start,format="%d-%b-%y",tz="EST")))/(3600*24)
		
		logLinCur=optim(c((log(repObs$Abund[1])-log(repObs$Abund[nrow(repObs)]))/(time[length(time)]-time[1]),log(max(repObs$Abund)),1),fitLogLinearDecay,N=log(repObs$Abund),time=time)
		logQuadCur=optim(c(((log(repObs$Abund[1])-log(repObs$Abund[nrow(repObs)]))/(time[length(time)]-time[1]))/10,(log(repObs$Abund[1])-log(repObs$Abund[nrow(repObs)]))/(time[length(time)]-time[1]),log(max(repObs$Abund)),1),fitLogQuadDecay,N=log(repObs$Abund),time=time)
		
		summ[counter,1]=strains[i]
		summ[counter,2]=reps[j]
		summ[counter,3:4]=logLinCur$par[1:2]
		summ[counter,5]=round(2*logLinCur$value+2*length(logLinCur$par),2)
		summ[counter,6:8]=logQuadCur$par[1:3]
		summ[counter,9]=round(2*logQuadCur$value+2*length(logQuadCur$par),2)
		
		summ[counter,10]=pchisq((2*logLinCur$value-2*logQuadCur$value),df=1,lower.tail=FALSE)
		
		counter=counter+1
		
		#plot(time,log(repObs$Abund),main=paste(strains[i]," rep ",reps[j]),ylim=c(0,20))
		#predTime=seq(0,max(time))
		#lines(predTime,logLinCur$par[2]-logLinCur$par[1]*predTime,lwd=2,lty=2)
		#lines(predTime,logQuadCur$par[3]-logQuadCur$par[1]*predTime^2-logQuadCur$par[2]*predTime,col='red',lwd=2,lty=2)
		}
	}
	
}

summ=summ[!is.na(summ[,1]),]
colnames(summ)=c('strain','rep','linearfit_K','linearfit_N0','linearfit_AIC','quadfit_m','quadfit_b','quadfit_N0','quadfit_AIC','LRT_pvalue')


#percent of tubes that are better fit by quadratic (AIC_linear > AIC_quadratic)
sum(as.numeric(summ[,5])>as.numeric(summ[,9]))/nrow(summ)*100	# 76%

# percent of tubes that are better based on likelihood ratio test
sum(as.numeric(summ[,10])<0.05)/nrow(summ)*100	# 63.7%

# add bonferroni correction
0.05/nrow(summ)	#  0.00044
sum(as.numeric(summ[,10])<0.0044)/nrow(summ)*100	# 40.7%

## histogram of likelihood ratio p-values
hist(log10(as.numeric(summ[,10])),breaks=seq(-20,0,1))
abline(v=log10(0.05),lwd=2,lty=2,col='red')
abline(v=log10(0.00044),lwd=2,lty=2,col='green')
legend('topleft',c('alpha=0.05','alpha=0.00044 (bonferroni)'),lty=2,col=c('red','green'),box.lty=0,lwd=2)

### confirm that likelihood ratio test is showing what deltaAIC between linear and quadratic would suggest
plot(log10(as.numeric(summ[,10])),as.numeric(summ[,9])-as.numeric(summ[,5]),xlab="log10(LRT pvalue)",ylab="deltaAIC (quadratic - linear)")
abline(v=log10(0.05),lwd=2,lty=2,col='red')
abline(v=log10(0.00044),lwd=2,lty=2,col='green')






plot(time,log(repObs$Abund))
predTime=seq(0,max(time))
lines(predTime,logLinCur$par[2]-logLinCur$par[1]*predTime,lwd=2,lty=2)
lines(predTime,logLinDecCur$par[3]-logLinDecCur$par[1]*predTime^2-logLinDecCur$par[2]*predTime,col='red',lwd=2,lty=2)

print(c("AIC of linear model:", round(2*logLinCur$value+2*length(logLinCur$par),2)))

print(c("AIC of quadratic model:",round(2*logQuadCur$value+2*length(logQuadCur$par))))

