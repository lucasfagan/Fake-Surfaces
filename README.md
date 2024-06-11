# Classification of Acyclic Fake Surfaces
This repository contains the data for the classification of acyclic cellular fake surfaces of complexity 1-4 and a partial classification of complexity 5. It also contains the  code used to generate the classification and check contractibility. 

## Fake Surface Data

## Classification Code
The code that generated the classification can be found in ``fakesurfaces_cla_6.py.`` In the filename, `cla` refers to the fact that it takes as input a command line argument corresponding to the 1-skeleton, and `6` that all 1-skeleta up to complexity 6 are written into the file.

## Checking Contractibility
For a given surface and maximal tree, we can put them at the top of this file as shown and then run ``sage check_contractibility.sage.`` It is also easy to modify the code to check multiple surfaces at once.  

## Generating One-Skeletons
We provide the code used to generate the one-skeleta in ``generate_one_skeleta.py``. Here we specify a complexity at the top, called ``SIZE_OF_MATRIX`` and then run the script.  
