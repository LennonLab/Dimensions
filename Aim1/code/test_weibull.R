rm(list = ls())
getwd()
setwd("~/GitHub/Dimensions/Aim1/")
## Load Data
obs <- read.csv("data/longtermdormancy_20151112_nocomments.csv", 
                header = TRUE, stringsAsFactors = FALSE)
## Adding 1 to deal with log(0) observations
obs$Abund <- as.numeric(obs$Colonies) * 10 ^ as.numeric(obs$Dilution) + 1
strains <- sort(unique(obs$Strain))
strains <- strains[table(obs$Strain)>10]
obs <- obs[obs$Strain%in%strains,]
summ <- matrix(NA,length(strains)*max(obs$Rep),7)

pdf('output/decayFitsWeibull.pdf') # Uncomment to create pdf that will plot data and fits
counter <- 1

for(i in 1:length(strains)){
  strainObs=obs[obs$Strain==strains[i],]
  
  reps=unique(strainObs$Rep)
  for(j in 1:length(reps)){
    
    repObs=strainObs[strainObs$Rep==reps[j],]
    
    if(nrow(repObs)>3){
      start=repObs[1,1]
      time=(as.numeric(strptime(repObs$Firstread_date,format="%d-%b-%y",tz="EST"))-
              as.numeric(strptime(start,format="%d-%b-%y",tz="EST")))/(3600*24)
      #logabund <- log10(repObs$Abund)
      time[time == 0] <- 0.1
      repObs["time"] <- time
      repObs["logabund"] <- log10(repObs$Abund)
      #print(repObs)

      # Initial parameters
      #A = 200 # Initial death (larger = slower) 
      #B = 1 # Bend (upper = 1 = first-order decay)
      #C = round(max(repObs$logabund),1) # intercept
      #Z = 6 # Error
      grids<-list(a=c(200),b=c(1),z=c(10))
      start<-list(a=NA,b=NA,c=round(max(repObs$logabund),1),z=NA)
      grid.starts<-as.matrix(expand.grid(grids))
      ncombos<-dim(grid.starts)[[1]]
      # cycle through each combo
      res.mat<-matrix(NA,nrow=ncombos,ncol=I(length(start)+1))
      res.mod<-list()
      for(k in 1:dim(grid.starts)[[1]]){
        #some how need to match grid parameters to start lists.
        mod.start<-as.list(grid.starts[k,])	
        new.start<-start
        new.start[names(start) %in% names(mod.start)]<-mod.start
        pscale<-as.numeric(new.start)
        names(pscale)<-names(new.start)
        fit <- mle2(minuslogl=logabund ~ dnorm(mean = c * (time / a)^(b-1) * exp(-1*(time/a)^b), sd = z), 
                                start = new.start, data = repObs, 
                                control=list(parscale=pscale, maxit=1000), 
                                method="Nelder-Mead", hessian = T)
                                #lower=c(a=0.0001, b=0.0001, c=2, z=0.0001), 
                                #upper=c(a=1000, b=1000, c=round(max(repObs$logabund),1), z=10),
                                #hessian = F)
                                
        #, method="L-BFGS-B", 
        #lower=c(a=0.0001, b=-10, c=2, z=0.0001)) 
        res.mat[k,]<-c(coef(fit),AIC(fit))		
        res.mod[[k]]<-fit
        print(res.mod)
      }
      colnames(res.mat)<-c(names(coef(fit)),"AIC")
      best.fit<-res.mod[[which(res.mat[,'AIC']==min(res.mat[,'AIC']))[1]]]
      summ[counter,1]=strains[i]
      summ[counter,2]=reps[j]
      # a
      summ[counter,3]=coef(best.fit)[1]
      # b
      summ[counter,4]=coef(best.fit)[2]
      # c
      summ[counter,5]=coef(best.fit)[3]
      # z
      summ[counter,6]=coef(best.fit)[4]
      summ[counter,7]=AIC(best.fit)
      
      #### add indicator of non-linearity (****) to plot title if quadratic model is better
      
      ### *** Comment/Uncomment following code to make pdf figs*** ###
      title=paste(strains[i]," rep ",reps[j])
      plot(repObs$time,log10(repObs$Abund),main=title,ylim=c(0,9))
      predTime=seq(0,max(repObs$time))
      lines(repObs$time, coef(best.fit)[3] * (repObs$time / coef(best.fit)[1])^(coef(best.fit)[2]-1) * exp(-1*(repObs$time/coef(best.fit)[1])^coef(best.fit)[2]), lwd=2,lty=2)
      #curve(coef(best.fit)[3] * (repObs$time / coef(best.fit)[1])^(coef(best.fit)[2]-1) * exp(-1*(repObs$time/coef(best.fit)[1])^coef(best.fit)[2]), 
      #      from = 0.1, to = 1000, add = TRUE, lty = 2, lwd = 4, col = "red") 
      # 			lines(predTime,logQuadCur$par[3]-logQuadCur$par[1]*predTime^2-logQuadCur$par[2]*predTime,col='red',lwd=2,lty=2)
      ### *** Comment/Uncomment above code to make pdf figs*** ###
      
      counter=counter+1
    }
  }
}

dev.off() 
summ=summ[!is.na(summ[,1]),]
colnames(summ)=c('strain','rep','a','b','c','z','AIC')
print(summ)
write.csv(summ,"data/weibull_results.csv")