from __future__ import division
import  matplotlib.pyplot as plt

import numpy as np
from random import choice
import scipy as sc
import scipy.spatial.distance as spd
import os
import sys


m1 = [1,1,1,1,1,1,1,1]
m2 = [1,1,1,1,1,1,1,1]
Dist = spd.braycurtis(m1, m2) # or jaccard
print 1-Dist


m1 = []
m2 = []

for i in range(10000):
    j = choice([0,1])
    m1.append(j)
    j = choice([0,1])
    m2.append(j)

Dist = spd.braycurtis(m1, m2) # or jaccard
print 1-Dist


m1 = [0]*28000
m2 = [0]*28000

for i in range(171100):
    j = choice(range(len(m1)))
    m1[j] += 1
    j = choice(range(len(m2)))
    m2[j] += 1

Dist = spd.braycurtis(m1, m2) # or jaccard
print 1-Dist

for i, val in enumerate(m1):
    if val > 1: m1[i] = 1
    if m2[i] > 1: m2[i] = 1


Dist = spd.braycurtis(m1, m2) # or jaccard
print 1-Dist
