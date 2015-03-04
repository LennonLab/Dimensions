# Notes from Sarah on files:

## XYponds.csv
PCNM analysis in R requires Cartesian coords. It was created by converting lat/long with geoXY in the SoDA package in R.

## HUC_10 and HUC_12
Differently scaled watersheds. HUC_10 is a large scale; HUC_12 is a larger scale. Kayla figured this out in GIS. The names of the real watersheds are in the second page. There is a map of these watersheds in the Indiana Ponds folder on the Lennon server.

## DNAtrans.csv 
Much of this had to do with getting rid of ponds for which I didn't have environmental data, singletons, and then doing a transformation.  See PNMC2.R for code, lines 39-49. # Select only OTUs with > 1 observations in >1 sample; Hellinger transformation; removed BC10, HNF221, YSF67.

## dorm_spatial.csv
This has to do specifically with variance partitioning between spatial and environmental matrices. It is a reduced matrix of only those spatial variables that are significantly correlated with the dormant matrix.  See Ponds_VarPartioning.R, line 190.  Similar file created for active matrix.

## INPonds.final.shared
Created by Mothur and contains # of reads for each of the OTUs created by Mothur at the 'unique', 97% similarity, and less than 97% similarity. All other data files descend from this file. Because of the way I set up the classification of OTUs on Mothur, the less than 97% similarity is not particularly useful.

## OTU.csv
Intermediary file made by Sarah Bray. It takes so long to load the shared file into memory that I (Sarah Bray) tried to minimize my load time by saving only the 97% similarity OTUs. My guess is that is what this file is from looking at the Excel file.  It looks like it contains all ponds for both DNA and cDNA for only the 97% similarity.