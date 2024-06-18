# Classification of Acyclic Fake Surfaces
This repository contains the data for the classification of acyclic cellular fake surfaces of complexity 1-4 and a partial classification of complexity 5: surfaces without small disks. It also contains the code used to generate the classification and check contractibility. See [our paper](https://arxiv.org/abs/2406.09439) (2406.09439) for more information.

## Fake Surface Data

All acyclic cellular fake surfaces of complexity 1-4 and complexity 5 without small disks can be found in ``fakesurfaces.csv.`` Each row represents an acyclic cellular fake surface. The first two columns represent the complexity followed the number of the 1-skeleton. Then, all the disk attaching maps are listed, where a negative entry denotes going along the edge in the opposite way. Following each disk, the two `Y/N` columns are the answers to whether the disk is embedded and has a trivial $T$-bundle, respectively. Information on how to construct vertex and edge labelings along with the ordering of the one-skeletons of a given complexity can be found in [our paper](https://arxiv.org/abs/2406.09439) or in the file `Surface_presentation_convention.pdf`.


## Classification Code
The code that generated the classification can be found in ``fakesurfaces_cla_6.py.`` In the filename, `cla` refers to the fact that it takes as input a command line argument corresponding to the 1-skeleton, and `6` that all 1-skeleta up to complexity 6 are written into the file.

## Checking Contractibility
For a given surface and maximal tree, we can put them at the top of this file as shown and then run ``sage check_contractibility.sage.`` It is also easy to modify the code to check multiple surfaces at once.  

## Generating One-Skeletons
We provide the code used to generate the one-skeleta in ``generate_one_skeleta.py``. Here we specify a complexity at the top, called ``SIZE_OF_MATRIX`` and then use the functions as named. Sample usage is provided at the bottom of the script.   
