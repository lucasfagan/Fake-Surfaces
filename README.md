# Classification of Acyclic Fake Surfaces
This repository contains the data for the classification of acyclic cellular fake surfaces of complexity 1-4 and a partial classification of complexity 5. It also contains the  code used to generate the classification and check contractibility. 

## Fake Surface Data

## Classification Code
The code that generated the classification can be found in ``fakesurfaces_cluster_cla_6.py.`` Cluster refers to the fact that it was prepared to run on the UCSB computing cluster, `cla` that it takes as input a command line argument corresponding to the 1-skeleton, and `6` that all 1-skeleta up to complexity 6 are written into the file.

## Checking Contractibility
For a given surface and maximal tree, we can put them at the top of this file as shown and then run ``sage check_contractibility.sage.`` It is also easy to modify the code to check multiple surfaces at once.  
