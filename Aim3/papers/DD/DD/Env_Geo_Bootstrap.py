from __future__ import division
import  matplotlib.pyplot as plt

import geopy
from geopy.distance import vincenty

import pandas as pd
import linecache
import numpy as np
import scipy as sc
import os
import sys

#import statsmodels.stats.api as sms
#import statsmodels.api as sm
import statsmodels.formula.api as smf
#from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table

mydir = os.path.expanduser("~/GitHub/Dimensions/Aim3/papers/DD")
mydir2 = os.path.expanduser("~/")


dat = pd.read_csv("~/GitHub/Dimensions/Aim3/DATA/EnvData/20130801_PondDataMod.csv", sep = ",", header = False)


dat2 = dat[dat['chla'] < 2000.0]
dat2 = dat2[dat2['pH'] > 1.0]
dat2 = dat2[dat2['Salinity'] > 0.0]
dat2 = dat2[dat2['TDS'] < 5.0]

results = pd.DataFrame()
trows = len(dat2.axes[0])
tcols = len(dat2.axes[1])

metrics = ['cityblock', 'euclidean', 'mahalanobis']

#SampleSize = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 47]
SampleSize = [4, 16, 12, 16]

for metric in metrics:
  RS =[]
  PVALS = []

  for size in SampleSize:
    rs = []
    pvals = []

    ct = 0
    while ct < 100:

      env = s.sample(n=size)

      lats = env['lat'] # latitudes (north and south)
      lons = env['long'] # longitudes (east and west)

      # Geographic Distances (Kilometers) Among Ponds
      latlong = as.matrix(cbind(env$long, env$lat))
      coord.dist = earth.dist(long.lat, dist = TRUE)
      #coord.dist = log(coord.dist)
      coord.dist[which(!is.finite(coord.dist))] = NA
      coord.dist.ls = liste(coord.dist, entry = "geo.dist")

      #env = subset(ponds, select=c(Depth:DON))
      #env.names = names(env)
      #env.dist = vegdist(env, "manhattan")

      salinity = as.numeric(env$"Salinity")
      spc = as.numeric(env$"SpC")
      tds = as.numeric(env$"TDS")
      doc = as.numeric(env$"DOC")
      env.dist = vegdist(cbind(salinity, spc, tds, doc), metric)

      # transform all distance matrices into list format:
      env.dist.ls = liste(env.dist, entry="env")

      df = data.frame(coord.dist.ls, env.dist.ls[,3])
      names(df)[3:4] = c("geo", "env")
      attach(df) #df = subset(df, struc != 0)

      r = cor(log10(df$env), df$geo)
      p = cor.test(log10(df$env), df$geo)$p.value

      rs = c(rs, r)
      pvals = c(pvals, p)
      ct = ct + 1
      }

    RS = c(RS, mean(rs))
    PVALS = c(PVALS, mean(pvals))
    }

  if(metric == 'manhattan'){
    results.df = cbind(RS, PVALS)
  }else{
    results.df = cbind(results.df, RS, PVALS)
  }
  }

cols = c("man_cor", "man_p", "euc_cor", "euc_p",
          "gow_cor", "gow_p", "mah_cor", "mah_p")
colnames(results.df) = cols
```


```{r, results = 'hide', message = FALSE, warning = FALSE}
plot.new()

results2.df = results.df
results2.df = na.omit(results2.df)
n = (SampleSize*(SampleSize-1))/2.0

par(mfrow=c(1, 2))

plot(n, results2.df[,2], xlab="sample size", ylab="p-val",
     main = "Significance vs. Sample size", col='SteelBlue',
     ylim = c(0, 0.7))

points(n, results2.df[,4], col='red4')
points(n, results2.df[,6], col='purple')
#points(n, results2.df[,8], col='green')

abline(0.05, 0.0, col='grey60')


plot(n, results2.df[,1], xlab="sample size", ylab="R",
     main = "R vs. Sample size", col='SteelBlue',
     ylim = c(0, 0.4))

points(n, results2.df[,3], col='red4')
points(n, results2.df[,5], col='purple')
#points(n, results2.df[,7], col='green')

abline(0.05, 0.0, col='grey60')
```
