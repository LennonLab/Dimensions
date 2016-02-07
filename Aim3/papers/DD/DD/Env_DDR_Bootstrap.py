from __future__ import division
import  matplotlib.pyplot as plt

import geopy
from geopy.distance import vincenty

import skbio
import skbio.diversity
import skbio.diversity.beta
from skbio.diversity import beta_diversity

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

EnvDat = pd.read_csv("~/GitHub/Dimensions/Aim3/DATA/EnvData/20130801_PondDataMod.csv", sep = ",", header = False)
Active = pd.read_csv("~/GitHub/Dimensions/Aim3/DATA/ForPython/ActiveComm.csv", sep = ",", header = False)
All = pd.read_csv("~/GitHub/Dimensions/Aim3/DATA/ForPython/AllComm.csv", sep = ",", header = False)
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
SampleSize = [4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45]
#SampleSize = [4, 8, 12, 20, 40]

BrayRS = []
JaccRS = []

BracPVALS = []
JaccPVALS = []

for size in SampleSize:
    print 'sample size:', size

    bray_rs = []
    jacc_rs = []

    bray_pvals = []
    jacc_pvals = []

    ct = 0
    while ct < 100:
        #inds = np.random.choice(range(1, trows+1), size, replace=False)

        GDists = []
        EnvDists = []
        BrayDists = []
        JaccDists = []

        for i in range(rows):
          row1 = env.iloc[[i]]
          lat1 = float(row1['lat']) # latitudes (north and south)
          long1 = float(row1['long']) # longitudes (east and west)
          env1 = row1.ix[:, 7:19]

          for j in range(rows):
            if j <= i: continue
            row2 = env.iloc[[j]]
            lat2 = float(row2['lat']) # latitudes (north and south)
            long2 = float(row2['long']) # longitudes (east and west)
            env2 = row2.ix[:, 7:19]

            # geographic distance
            geo_dist = vincenty((lat1, long1), (lat2, long2)).km
            GDists.append(geo_dist)

            # environmental distances
            eu_dist = spd.euclidean(env1, env2) # Euclidean distance
            EuDists.append(eu_dist)

            man_dist = spd.cityblock(env1, env2) # Manhattan distance
            ManDists.append(man_dist)

            ham_dist = spd.hamming(env1, env2) # Square Euclidean distance
            HamDists.append(ham_dist)

            squ_dist = spd.sqeuclidean(env1, env2) # Square Euclidean distance
            SquDists.append(squ_dist)

            cos_dist = spd.cosine(env1, env2) # cosine distance
            CosDists.append(cos_dist)

            cor_dist = spd.correlation(env1, env2) # correlation distance
            CorDists.append(cor_dist)

        # get correlation coefficient and p-value
        r, p = sc.stats.pearsonr(GDists, EuDists)
        euc_rs.append(r)
        euc_pvals.append(p)

        r, p = sc.stats.pearsonr(GDists, SquDists)
        squ_rs.append(r)
        squ_pvals.append(p)

        r, p = sc.stats.pearsonr(GDists, ManDists)
        man_rs.append(r)
        man_pvals.append(p)

        r, p = sc.stats.pearsonr(GDists, CorDists)
        cor_rs.append(r)
        cor_pvals.append(p)

        r, p = sc.stats.pearsonr(GDists, HamDists)
        ham_rs.append(r)
        ham_pvals.append(p)

        r, p = sc.stats.pearsonr(GDists, CosDists)
        cos_rs.append(r)
        cos_pvals.append(p)

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

    SquRS.append(np.mean(squ_rs))
    SquPVALS.append(np.mean(squ_pvals))

    CosRS.append(np.mean(cos_rs))
    CosPVALS.append(np.mean(cos_pvals))


print "generating figure"
fig = plt.figure()
fig.add_subplot(2, 2, 2)

SampleSize = np.array(SampleSize)
SampleSize = np.log10((SampleSize*(SampleSize - 1))/2)

plt.plot(SampleSize, SquPVALS, color = '0.2', alpha= 0.6 , linewidth = 2, label='Square Euc.')
plt.plot(SampleSize, ManPVALS, color = 'SteelBlue', alpha= 0.6 , linewidth=2, label='Manhattan')
plt.plot(SampleSize, EucPVALS, color = 'm', alpha= 0.6 , linewidth=2, label='Euclidean')
plt.plot(SampleSize, HamPVALS, color = '0.7', alpha= 0.6 , linewidth=2, label='Hamming')
plt.plot(SampleSize, CorPVALS, color = 'Limegreen', alpha= 0.6 , linewidth=2, label='Correlation')
plt.plot(SampleSize, CosPVALS, color = 'red', alpha= 0.6 , linewidth=2, label='Cosine')

plt.plot([min(SampleSize), max(SampleSize)], [0.05, 0.05], c='0.2', ls='--')
plt.legend(bbox_to_anchor=(-0.03, 1.07, 2.46, .3), loc=10, ncol=3, mode="expand",prop={'size':12})

plt.ylim(0, 0.2)
plt.ylabel('$p$'+'-value', fontsize=14)
plt.xlabel('sample size', fontsize=14)

fig.add_subplot(2, 2, 2)
plt.plot(SampleSize, SquRS, color = '0.2', alpha= 0.6 , linewidth = 2)
plt.plot(SampleSize, ManRS, color = 'SteelBlue', alpha= 0.6 ,  linewidth=2)
plt.plot(SampleSize, EucRS, color = 'm', alpha= 0.6 , linewidth=2)
plt.plot(SampleSize, HamRS, color = '0.7', alpha= 0.6 , linewidth=2)
plt.plot(SampleSize, CorRS, color = 'Limegreen', alpha= 0.6 , linewidth=2)
plt.plot(SampleSize, CosRS, color = 'red', alpha= 0.6 , linewidth=2)

plt.ylabel('Correlation coefficient', fontsize=12)
plt.xlabel('sample size', fontsize=12)

plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir+'/DD/figs/SampleSizeDependency.png', dpi=300, bbox_inches = "tight")
#plt.show()
