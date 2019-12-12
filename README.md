# Algebraic Characterization of Essential Matrices and Their Averaging in Multiview Settings (ICCV 2019)

Yoni Kasten* Amnon Geifman* Meirav Galun Ronen Basri (*equal contribution)


MATLAB Mex Version 1.0 (2019-12-12)
Copyright (C) Yoni Kasten and Amnon Geifman, Weizmann Institute, 2019.
Licensed for noncommercial research use only.


## Background

The code retrieves camera matrices given a set of pairwise essential matrices.

For more information see:

[[Paper]](http://openaccess.thecvf.com/content_ICCV_2019/papers/Kasten_Algebraic_Characterization_of_Essential_Matrices_and_Their_Averaging_in_Multiview_ICCV_2019_paper.pdf)
Please cite this paper if you use this code in an academic publication.
```
@InProceedings{Kasten_2019_ICCV,
	author = {Kasten, Yoni and Geifman, Amnon and Galun, Meirav and Basri, Ronen},
	title = {Algebraic Characterization of Essential Matrices and Their Averaging in Multiview Settings},
	booktitle = {The IEEE International Conference on Computer Vision (ICCV)},
	month = {October},
	year = {2019}
}
```

## Installation


mex/C++ code:
In order to use the code it is necessary to compile the mex functions.
We supply compiled versions for Windows 

Compiling commands:
```
mex formexgraph6.c
mex buildTripletGraph.c
```



## Use

We supply end-to-end code for the pipeline described in the paper.
Important function:
```
runAll.m - Script to run the pipline on the specified data sets. Outputs a table with the results and a .mat file that contains the camera parameters. The script generates table.csv.
GetTableFiltered.m - Script to get the results on the reconstructed cameras as apears in the paper. The BA of Theia filters some of the cameras before BA. We keep the list at Final_names. The script uses the results which saved into 'results' folder. (need to run "runAll.m" before). The script generates tableFinal.csv.



```
Important variables:

```
pointMatchesInliers1- a matrix that contains the number of inliers between each pair of frames
EN- the multi-view essential matrix (measurements matrix).
Hmat- a block matrix where the i,j 3x3 block is the relative rotation between the i'th and the j'th frames.
TijMat- a cell array where the i,j entry is the relative translation between the i'th and the j'th frames.
R_gt,T_gt,K_gt - the estimated ground truth rotation translation and calibration (respectively)
M- tracks matrix 

```

## Acknowledgement 
We use the following 3rdparties code:
```
SimpleTransScaleRemove.m
GlobalSOdCorrectLeft
getBestError
```
This files are from LUD's code. (Author: Onur Ozyesil). 
We use them to find the best alignment of the results with the ground truth.


## Contact 
For any query, contact : 
Yoni Kasten, Amnon Geifman 
Weizmann Institute of Science
{yoni.kasten,amnon.geifman}@weizmann.ac.il

## License
   This software is provided under the provisions of the Lesser GNU Public License (LGPL). 
   see: http://www.gnu.org/copyleft/lesser.html.

   This software can be used only for research purposes, you should cite
   the aforementioned papers in any resulting publication.

   The Software is provided "as is", without warranty of any kind.




## Version History


* Version 1.0 (2019-12-12)
   Initial Release
