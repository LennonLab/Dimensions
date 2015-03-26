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
fitLogExpDecay<-function(p,N,time){
	k=p[1]
	N0=p[2]
	sd=exp(p[3])
	
	
	Nt=N0-k*time
	
	-sum(dnorm(N,Nt,sd,log=TRUE))
}

# likelihood function for quadratic fit to log transformed data
fitLogDecDecay<-function(p,N,time){
	m=p[1]
	b=p[2]
	N0=p[3]
	sd=exp(p[4])

	Nt=N0-m*time^2-b*time
	
	-sum(dnorm(N,Nt,sd,log=TRUE))
}


strains=sort(unique(obs$Strain))
#only fit strains with more than 10 observations
strains=strains[table(obs$Strain)>10]
obs=obs[obs$Strain%in%strains,]

#pdf('decayFits.pdf')

summ=matrix(NA,length(strains)*max(obs$Rep),9)
counter=1
for(i in 1:length(strains)){
	strainObs=obs[obs$Strain==strains[i],]
	
	reps=unique(strainObs$Rep)
	for(j in 1:length(reps)){
		
		repObs=strainObs[strainObs$Rep==reps[j],]

		if(nrow(repObs)>3){
			
		start=repObs[1,1]
		
		time=(as.numeric(strptime(repObs$Firstread_date,format="%d-%b-%y",tz="EST"))-as.numeric(strptime(start,format="%d-%b-%y",tz="EST")))/(3600*24)
		
		logLinCur=optim(c((log(repObs$Abund[1])-log(repObs$Abund[nrow(repObs)]))/(time[length(time)]-time[1]),log(max(repObs$Abund)),1),fitLogExpDecay,N=log(repObs$Abund),time=time)
		logLinDecCur=optim(c(((log(repObs$Abund[1])-log(repObs$Abund[nrow(repObs)]))/(time[length(time)]-time[1]))/10,(log(repObs$Abund[1])-log(repObs$Abund[nrow(repObs)]))/(time[length(time)]-time[1]),log(max(repObs$Abund)),1),fitLogDecDecay,N=log(repObs$Abund),time=time)
		
		summ[counter,1]=strains[i]
		summ[counter,2]=reps[j]
		summ[counter,3:4]=logLinCur$par[1:2]
		summ[counter,5]=round(2*logLinCur$value+2*length(logLinCur$par),2)
		summ[counter,6:8]=logLinDecCur$par[1:3]
		summ[counter,9]=round(2*logLinDecCur$value+2*length(logLinDecCur$par))
		
		counter=counter+1
		
		#plot(time,log(repObs$Abund),main=paste(strains[i]," rep ",reps[j]),ylim=c(0,20))
		predTime=seq(0,max(time))
		#lines(predTime,logLinCur$par[2]-logLinCur$par[1]*predTime,lwd=2,lty=2)
		#lines(predTime,logLinDecCur$par[3]-logLinDecCur$par[1]*predTime^2-logLinDecCur$par[2]*predTime,col='red',lwd=2,lty=2)

		
		}
	}
	
}

summ=summ[!is.na(summ[,1]),]
colnames(summ)=c('strain','rep','linearfit_K','linearfit_N0','linearfit_AIC','quadfit_m','quadfit_b','quadfit_N0','quadfit_AIC')


#percent of tubes that are better fit by quadratic
sum(as.numeric(summ[,5])>as.numeric(summ[,9]))/nrow(summ)*100


plot(time,log(repObs$Abund))
predTime=seq(0,max(time))
lines(predTime,logLinCur$par[2]-logLinCur$par[1]*predTime,lwd=2,lty=2)
lines(predTime,logLinDecCur$par[3]-logLinDecCur$par[1]*predTime^2-logLinDecCur$par[2]*predTime,col='red',lwd=2,lty=2)

print(c("AIC of linear model:", round(2*logLinCur$value+2*length(logLinCur$par),2)))

print(c("AIC of quadratic model:",round(2*logLinDecCur$value+2*length(logLinDecCur$par))))

