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

import skbio
import skbio.stats.ordination as ord

import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table


mydir = os.path.expanduser("~/GitHub/Dimensions/Aim3/papers/DD")
mydir2 = os.path.expanduser("~/")

EnvDat = pd.read_csv("~/GitHub/Dimensions/Aim3/DATA/EnvData/20130801_PondDataMod.csv", sep = ",", header = 0)
Active = pd.read_csv("~/GitHub/Dimensions/Aim3/DATA/ForPython/ActiveComm.csv", sep = ",", header = 0)
All = pd.read_csv("~/GitHub/Dimensions/Aim3/DATA/ForPython/AllComm.csv", sep = ",", header = 0)

dat2 = EnvDat[EnvDat['chla'] < 2000.0]
dat2 = dat2[dat2['pH'] > 1.0]
dat2 = dat2[dat2['Salinity'] > 0.0]
dat2 = dat2[dat2['TDS'] < 5.0]

trows = len(dat2.axes[0])
tcols = len(dat2.axes[1])
ColNames = list(dat2.columns.values)

for i, column in enumerate(dat2):
    if i >= 6:
        #print ColNames[i], i
        #print dat2.ix[:,i]
        dat2.ix[:,i] = (dat2.ix[:,i] - np.mean(dat2.ix[:,i]))/np.std(dat2.ix[:,i])


GeoDists = []
EnvDists = []
ActDists = []
AllDists = []


for i in range(trows):
    row1 = dat2.iloc[[i]]
    lat1 = float(row1['lat']) # latitudes (north and south)
    long1 = float(row1['long']) # longitudes (east and west)

    env1 = row1.ix[:, 7:19]
    act_row1 = Active.iloc[[i]]
    all_row1 = All.iloc[[i]]

    for j in range(trows):
        if j <= i: continue
        row2 = dat2.iloc[[j]]
        lat2 = float(row2['lat']) # latitudes (north and south)
        long2 = float(row2['long']) # longitudes (east and west)

        env2 = row2.ix[:, 7:19]
        act_row2 = Active.iloc[[j]]
        all_row2 = All.iloc[[j]]

        # geographic distance
        geo_dist = vincenty((lat1, long1), (lat2, long2)).km
        GeoDists.append(geo_dist)

        # environmental distance
        env_dist = spd.correlation(env1, env2) # correlation distance
        EnvDists.append(env_dist)

        # beta-diversity
        act_com_dist = spd.braycurtis(act_row1, act_row2) # or jaccard
        ActDists.append(1 - act_com_dist)

        all_com_dist = spd.braycurtis(all_row1, all_row2)
        AllDists.append(1 - all_com_dist)


fig = plt.figure()

ax = fig.add_subplot(2, 2, 1)
m, b, r, p, stder = sc.stats.linregress(GeoDists, ActDists)
plt.scatter(GeoDists, ActDists, color = '0.4', alpha= 0.6, s = 5, linewidth = 2, label='Active, Geographic')

plt.ylim(0, 1.0)
ax.text(0, 0.8, 'Active, Geographic', fontsize=12, color = '0.4')
plt.ylabel('Bray-Curtis', fontsize=14)
plt.xlabel('Kilometers', fontsize=14)


ax = fig.add_subplot(2, 2, 2)
m, b, r, p, stder = sc.stats.linregress(EnvDists, ActDists)
plt.scatter(EnvDists, ActDists, color = 'SteelBlue', alpha= 0.6, s=5, linewidth=2, label='Active, Environmental')

plt.ylim(0.0, 1.0)
ax.text(0, 0.8, 'Active, Environmental', fontsize=12, color = 'SteelBlue')
plt.ylabel('Bray-Curtis', fontsize=14)
plt.xlabel('CCA1', fontsize=14)


ax = fig.add_subplot(2, 2, 3)
m, b, r, p, stder = sc.stats.linregress(GeoDists, AllDists)
plt.scatter(GeoDists, AllDists, color = 'm', alpha= 0.6, linewidth=2, s=5, label='All, Geographic')

plt.ylim(0, 1)
ax.text(0, 0.8, 'All, Geographic', fontsize=12, color = 'm')
plt.ylabel('Bray-Curtis', fontsize=14)
plt.xlabel('Kilometers', fontsize=14)


ax = fig.add_subplot(2, 2, 4)
m, b, r, p, stder = sc.stats.linregress(EnvDists, AllDists)
plt.scatter(EnvDists, AllDists, color = 'DarkTurquoise', alpha= 0.2 , s=5, linewidth=2, label='All, Environmental')

plt.ylim(0, 1)
ax.text(0, 0.8, 'All, Environmental', fontsize=12, color = 'DarkTurquoise')
plt.ylabel('Bray-Curtis', fontsize=14)
plt.xlabel('CCA1', fontsize=14)


plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir+'/DD/figs/DDRs.png', dpi=300, bbox_inches = "tight")
#plt.show()
