from __future__ import division
import  matplotlib.pyplot as plt

import geopy
from geopy.distance import vincenty

import pandas as pd
import linecache
import numpy as np
import scipy as sc
import scipy.spatial.distance as spd
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

ColNames = list(dat2.columns.values)

for i, column in enumerate(dat2):
    if i >= 6:
        #print ColNames[i], i
        #print dat2.ix[:,i]
        dat2.ix[:,i] = (dat2.ix[:,i] - np.mean(dat2.ix[:,i]))/np.std(dat2.ix[:,i])


metrics = ['Euclidean', 'Hamming']

#SampleSize = [4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45]
SampleSize = [4, 8, 12, 20, 40, 47]

ManRS = []
HamRS = []
CorRS =[]
EucRS =[]

ManPVALS = []
HamPVALS = []
CorPVALS = []
EucPVALS = []


for size in SampleSize:
    print 'sample size:', size
    man_rs = []
    ham_rs = []
    cor_rs = []
    euc_rs = []

    man_pvals = []
    ham_pvals = []
    cor_pvals = []
    euc_pvals = []

    ct = 0
    while ct < 100:
        GDists = []
        ManDists = []
        HamDists = []
        EuDists = []
        CorDists = []

        env = dat2.sample(n=size)
        rows = len(env.axes[0])

        for i in range(rows):
          row1 = env.iloc[[i]]
          lat1 = float(row1['lat']) # latitudes (north and south)
          long1 = float(row1['long']) # longitudes (east and west)
          env1 = row1.ix[:, 5:19]

          for j in range(rows):
            if j <= i: continue
            row2 = env.iloc[[j]]
            lat2 = float(row2['lat']) # latitudes (north and south)
            long2 = float(row2['long']) # longitudes (east and west)
            env2 = row2.ix[:, 5:19]

            # geographic distance
            geo_dist = vincenty((lat1, long1), (lat2, long2)).km
            GDists.append(geo_dist)

            # environmental distances
            eu_dist = spd.euclidean(env1, env2) # Euclidean distance
            EuDists.append(eu_dist)

            man_dist = spd.sqeuclidean(env1, env2) # Manhattan distance
            ManDists.append(man_dist)

            ham_dist = spd.hamming(env1, env2) # hamming distance
            HamDists.append(ham_dist)

            cor_dist = spd.correlation(env1, env2) # correlation distance
            CorDists.append(cor_dist)

        # get correlation coefficient and p-value
        r1, p1 = sc.stats.pearsonr(GDists, EuDists)
        r2, p2 = sc.stats.pearsonr(GDists, ManDists)
        r3, p3 = sc.stats.pearsonr(GDists, CorDists)
        r4, p4 = sc.stats.pearsonr(GDists, HamDists)

        man_rs.append(r1)
        man_pvals.append(p1)
        cor_rs.append(r2)
        cor_pvals.append(p2)
        euc_rs.append(r3)
        euc_pvals.append(p3)
        ham_rs.append(r4)
        ham_pvals.append(p4)

        ct += 1

    # get average cc and p-val
    ManRS.append(np.mean(man_rs))
    ManPVALS.append(np.mean(man_pvals))

    HamRS.append(np.mean(ham_rs))
    HamPVALS.append(np.mean(ham_pvals))

    CorRS.append(np.mean(cor_rs))
    CorPVALS.append(np.mean(cor_pvals))

    EucRS.append(np.mean(euc_rs))
    EucPVALS.append(np.mean(euc_pvals))


print "generating figure"
fig = plt.figure()
fig.add_subplot(2, 2, 1)

SampleSize = np.array(SampleSize)
SampleSize = np.log10((SampleSize*(SampleSize - 1))/2)

size = 60

plt.scatter(SampleSize, ManPVALS, color = '0.4', alpha= 0.6 , s = size, linewidths=2, edgecolor='0.2', label='Manhattan')
plt.scatter(SampleSize, CorPVALS, color = 'cyan', alpha= 0.6 , s = size, linewidths=2, edgecolor='Steelblue', label='Correlation')
plt.scatter(SampleSize, EucPVALS, color = 'm', alpha= 0.6 , s = size, linewidths=2, edgecolor='m', label='Euclidean')
plt.scatter(SampleSize, HamPVALS, color = '0.7', alpha= 0.6 , s = size, linewidths=2, edgecolor='0.5', label='Hamming')
plt.legend(bbox_to_anchor=(-0.03, 1.07, 2.46, .3), loc=10, ncol=2, mode="expand",prop={'size':12})

plt.ylabel('$p$'+'-value', fontsize=14)
plt.xlabel('sample size', fontsize=14)

fig.add_subplot(2, 2, 2)
plt.scatter(SampleSize, ManRS, color = '0.4', alpha= 0.6 , s = size, linewidths=2, edgecolor='0.2', label='Manhattan')
plt.scatter(SampleSize, CorRS, color = 'cyan', alpha= 0.6 , s = size, linewidths=2, edgecolor='Steelblue', label='Mahalanobis')
plt.scatter(SampleSize, EucRS, color = 'm', alpha= 0.6 , s = size, linewidths=2, edgecolor='m', label='Euclidean')
plt.scatter(SampleSize, HamRS, color = '0.7', alpha= 0.6 , s = size, linewidths=2, edgecolor='0.5', label='Hamming')

plt.ylabel('Correlation coefficient', fontsize=12)
plt.xlabel('sample size', fontsize=12)

plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir+'/DD/figs/SampleSizeDependency.png', dpi=300, bbox_inches = "tight")
#plt.show()
