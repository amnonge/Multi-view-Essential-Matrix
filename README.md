# Algebraic Characterization of Essential Matrices and Their Averaging in Multiview Settings (ICCV 2019)

Yoni Kasten* Amnon Geifman* Meirav Galun Ronen Basri (*equal contribution)


MATLAB Mex Version 1.0 (2019-06-06)
Copyright (C) Yoni Kasten and Amnon Geifman, Weizmann Institute, 2019.
Licensed for noncommercial research use only.


## Background

The code retrieves camera matrices given a set of pairwise essential matrices.

For more information see:

[[arXiv]](https://arxiv.org/abs/1904.02663)
Please cite these paper if you use this code in an academic publication.
```
@article{kasten2019algebraic,
  title={Algebraic Characterization of Essential Matrices and Their Averaging in Multiview Settings},
  author={Kasten, Yoni and Geifman, Amnon and Galun, Meirav and Basri, Ronen},
  journal={arXiv preprint arXiv:1904.02663},
  year={2019}
}
```

## Installation


mex/C++ code:
In order to use the code it is necessary to compile the mex functions.
We supply compiled versions for Windows and Linux.





## Use

We supply end-to-end code for the pipeline described in the paper.
Important function:
```
runall.m- Script to run the pipline on the specified data sets. Outputs a table with the results and a .mat file that contains the camera parameters

GetTableAfterBA.m - Script to get the results on the reconstructed cameras as apears in the paper. This subset of camera are the reconstructed cameras got from the BA. The script uses the results which saved into 'results' folder.



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
3rdparty\fromPPSFM - evaluation and self calibration code from "Practical projective structure
from motion (p2sfm)" (ICCV 2017). Their self calibration code implements the paper "Autocalibration via rank-constrained estimation of the
absolute quadric" (CVPR 2007), and uses the libraries "GloptiPoly 3" and "SeDuMi 1.3"

3rdparty\vgg_code - implementation of basic functions from the book "Multiple View Geometry in Computer Vision" (2004) downloaded from:
https://www.robots.ox.ac.uk/~vgg/hzbook/code/

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


* Version 1.0 (2019-06-06)
   Initial Release
