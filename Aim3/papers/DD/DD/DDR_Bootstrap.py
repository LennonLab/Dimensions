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

EnvDat = pd.read_csv("~/GitHub/Dimensions/Aim3/DATA/EnvData/20130801_PondDataMod.csv", sep = ",", header = 0)
Active = pd.read_csv("~/GitHub/Dimensions/Aim3/DATA/ForPython/ActiveComm.csv", sep = ",", header = 0)
All = pd.read_csv("~/GitHub/Dimensions/Aim3/DATA/ForPython/AllComm.csv", sep = ",", header = 0)
#ColNames = list(Active.columns.values)
#print ColNames

dat2 = EnvDat[EnvDat['chla'] < 2000.0]
dat2 = dat2[dat2['pH'] > 1.0]
dat2 = dat2[dat2['Salinity'] > 0.0]
dat2 = dat2[dat2['TDS'] < 5.0]


results = pd.DataFrame()
trows = len(dat2.axes[0])
tcols = len(dat2.axes[1])
allrows = len(All.axes[0])
actrows = len(Active.axes[0])
if trows != allrows or trows != actrows: sys.exit()

ColNames = list(dat2.columns.values)

for i, column in enumerate(dat2):
    if i >= 6:
        #print ColNames[i], i
        #print dat2.ix[:,i]
        dat2.ix[:,i] = (dat2.ix[:,i] - np.mean(dat2.ix[:,i]))/np.std(dat2.ix[:,i])

#sys.exit()
#SampleSize = [4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45]
SampleSize = [4, 8, 10, 20, 30]

ActEnvRS = []
ActGeoRS = []
AllEnvRS = []
AllGeoRS = []

ActEnvPVALS = []
ActGeoPVALS = []
AllEnvPVALS = []
AllGeoPVALS = []

ActEnvSlope = []
ActGeoSlope = []
AllEnvSlope = []
AllGeoSlope = []


for size in SampleSize:
    print 'sample size:', size

    _ActEnvRS = []
    _ActGeoRS = []
    _AllEnvRS = []
    _AllGeoRS = []

    _ActEnvPVALS = []
    _ActGeoPVALS = []
    _AllEnvPVALS = []
    _AllGeoPVALS = []

    _ActEnvSlope = []
    _ActGeoSlope = []
    _AllEnvSlope = []
    _AllGeoSlope = []

    ct = 0
    while ct < 10:
        inds = np.random.choice(range(0, trows), size, replace=False)

        GeoDists = []
        EnvDists = []
        AllDists = []
        ActDists = []

        for i, val1 in enumerate(inds):
          row1 = dat2.iloc[[val1]]
          lat1 = float(row1['lat']) # latitudes (north and south)
          long1 = float(row1['long']) # longitudes (east and west)

          env1 = row1.ix[:, 7:19]
          act_row1 = Active.iloc[[val1]]
          all_row1 = All.iloc[[val1]]

          for j, val2 in enumerate(inds):
            if j <= i: continue
            row2 = dat2.iloc[[val2]]
            lat2 = float(row2['lat']) # latitudes (north and south)
            long2 = float(row2['long']) # longitudes (east and west)

            env2 = row2.ix[:, 7:19]
            act_row2 = Active.iloc[[val2]]
            all_row2 = All.iloc[[val2]]

            # geographic distance
            geo_dist = vincenty((lat1, long1), (lat2, long2)).km
            GeoDists.append(geo_dist)

            # environmental distance
            env_dist = spd.correlation(env1, env2) # correlation distance
            EnvDists.append(env_dist)

            # beta-diversity
            act_com_dist = spd.braycurtis(act_row1, act_row2) # or jaccard
            ActDists.append(act_com_dist)

            all_com_dist = spd.braycurtis(all_row1, all_row2)
            AllDists.append(all_com_dist)


        # get correlation coefficient and p-value
        m, b, r, p, stder = sc.stats.linregress(GeoDists, ActDists)
        _ActGeoRS.append(r)
        _ActGeoPVALS.append(p)
        _ActGeoSlope.append(m)
        _ActGeoINT.append(b)


        m, b, r, p, stder = sc.stats.linregress(EnvDists, ActDists)
        _ActEnvRS.append(r)
        _ActEnvPVALS.append(p)
        _ActEnvSlope.append(m)
        _ActEnvINT.append(b)

        m, b, r, p, stder = sc.stats.linregress(GeoDists, AllDists)
        _AllGeoRS.append(r)
        _AllGeoPVALS.append(p)
        _AllGeoSlope.append(m)
        _AllGeoINT.append(b)

        m, b, r, p, stder = sc.stats.linregress(EnvDists, AllDists)
        _AllEnvRS.append(r)
        _AllEnvPVALS.append(p)
        _AllEnvSlope.append(m)
        _AllEnvINT.append(b)

        ct += 1

    # get average cc and p-val
    ActGeoRS.append(np.mean(_ActGeoRS))
    ActGeoPVALS.append(np.mean(_ActGeoPVALS))
    ActGeoSlope.append(np.mean(_ActGeoSlope))
    ActGeoINT.append(np.mean(_ActGeoINT))

    ActEnvRS.append(np.mean(_ActEnvRS))
    ActEnvPVALS.append(np.mean(_ActEnvPVALS))
    ActEnvSlope.append(np.mean(_ActEnvSlope))
    ActEnvINT.append(np.mean(_ActEnvINT))

    AllGeoRS.append(np.mean(_AllGeoRS))
    AllGeoPVALS.append(np.mean(_AllGeoPVALS))
    AllGeoSlope.append(np.mean(_AllGeoSlope))
    AllGeoINT.append(np.mean(_AllGeoINT))

    AllEnvRS.append(np.mean(_AllEnvRS))
    AllEnvPVALS.append(np.mean(_AllEnvPVALS))
    AllEnvSlope.append(np.mean(_AllEnvSlope))
    AllEnvINT.append(np.mean(_AllEnvINT))


print "generating figure"
fig = plt.figure()
fig.add_subplot(2, 2, 1)

SampleSize = np.array(SampleSize)
SampleSize = np.log10((SampleSize*(SampleSize - 1))/2)

plt.plot(SampleSize, ActGeoPVALS, color = '0.2', alpha= 0.6 , linewidth = 2, label='Active Geo')
plt.plot(SampleSize, ActEnvPVALS, color = 'SteelBlue', alpha= 0.6 , linewidth=2, label='Active Env')
plt.plot(SampleSize, AllGeoPVALS, color = 'm', alpha= 0.6 , linewidth=2, label='All Geo')
plt.plot(SampleSize, AllEnvPVALS, color = '0.7', alpha= 0.6 , linewidth=2, label='All Env')

plt.plot([min(SampleSize), max(SampleSize)], [0.05, 0.05], c='0.2', ls='--')
plt.legend(bbox_to_anchor=(-0.03, 1.07, 2.46, .3), loc=10, ncol=2, mode="expand",prop={'size':12})

plt.ylim(0, 0.2)
plt.ylabel('$p$'+'-value', fontsize=14)
plt.xlabel('sample size', fontsize=14)

fig.add_subplot(2, 2, 2)
plt.plot(SampleSize, ActGeoRS, color = '0.2', alpha= 0.6 , linewidth = 2)
plt.plot(SampleSize, ActEnvRS, color = 'SteelBlue', alpha= 0.6 , linewidth=2)
plt.plot(SampleSize, AllGeoRS, color = 'm', alpha= 0.6 , linewidth=2)
plt.plot(SampleSize, AllEnvRS, color = '0.7', alpha= 0.6 , linewidth=2)

plt.ylabel('Correlation coefficient', fontsize=12)
plt.xlabel('sample size', fontsize=12)


fig.add_subplot(2, 2, 3)

plt.plot(SampleSize, ActGeoSlope, color = '0.2', alpha= 0.6 , linewidth = 2)
plt.plot(SampleSize, ActEnvSlope, color = 'SteelBlue', alpha= 0.6 , linewidth=2)
plt.plot(SampleSize, AllGeoSlope, color = 'm', alpha= 0.6 , linewidth=2)
plt.plot(SampleSize, AllEnvSlope, color = '0.7', alpha= 0.6 , linewidth=2)

plt.ylabel('DDR slope', fontsize=14)
plt.xlabel('sample size', fontsize=14)

fig.add_subplot(2, 2, 4)

plt.plot(SampleSize, ActGeoSlope, color = '0.2', alpha= 0.6 , linewidth = 2)
plt.plot(SampleSize, ActEnvSlope, color = 'SteelBlue', alpha= 0.6 , linewidth=2)
plt.plot(SampleSize, AllGeoSlope, color = 'm', alpha= 0.6 , linewidth=2)
plt.plot(SampleSize, AllEnvSlope, color = '0.7', alpha= 0.6 , linewidth=2)

plt.ylabel('DDR slope', fontsize=14)
plt.xlabel('sample size', fontsize=14)

plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir+'/DD/figs/DDR-SampleSizeDependency.png', dpi=300, bbox_inches = "tight")
#plt.show()
