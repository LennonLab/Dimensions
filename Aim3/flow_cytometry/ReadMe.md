# flow_cytometry

This directory contains data files produced by the flow analysis of each IN Pond sample

## Contents
### Directories and files

**INPonds_flowdat.csv**: summary output data from each IN pond sample.

Column descriptions

 * *pond*: IN pond ID
 * *sample*: flow sample ID
 * *meanDNA*: mean DNA fluorescence
 * *meanRNA: mean RNA fluorescence
 * *mean.RNADNA*: mean RNA/DNA ratio
 * *sd.RNADNA* : standard deviation of RNA/DNA ratio
 * *median.RNADNA* median RNA/DNA ratio
 * *N.totdens*: adjusted total density (cells/mL) of bacteria in community
 * *N.deaddens*: adjusted density (cells/mL) of dead bacteria in community
 * *N.livedens*: adjusted density (cells/mL) of live bacteria in community


**event_counts/**: This directory contains the RNA/DNA ratio as well as the RNA and DNA fluorescence per event within each pond. Files are named according to sample name as described above.

**figures/**: This directory contains .png images of the RNA/DNA histograms


