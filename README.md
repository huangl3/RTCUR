# Mode-wise Tensor Decomposition

This repository is for our paper:

[1] HanQin Cai, Zehan Chao, Longxiu Huang, and Deanna Needell. <a href=https://doi.org/10.1137/23M1574282> Robust Tensor CUR: Rapid Low-Tucker-Rank Tensor Recovery with Sparse Corruptions</a>. *SIAM Journal on Imaging Sciences*, 17(1): 225–247, 2024. 

###### To display math symbols properly, one may have to install a MathJax plugin. For example, [MathJax Plugin for Github](https://chrome.google.com/webstore/detail/mathjax-plugin-for-github/ioemnmodlmafdkllaclgeombjnmnbima?hl=en).


## Introduction
A fast non-convex optimization algorithm, called RTCUR, for tensor robust principal component analysis (TRPCA) problem. We provide four verisons of this generalization.


## Environment
This repo is developed with <a href=https://gitlab.com/tensors/tensor_toolbox/-/releases/v3.1>Tensor Toolbox v3.1</a>. A future verison of this toolbox may also increase the peformance of our code; however, we cannot guarantee their compatibility.


## Syntex
All functions in the algorithm folder follow the same syntax. For example:
#### RTCUR with fixed Chidori CUR
```
[L_core,X_sub_mat,timer,err ] = RTCUR_fc(D, r, para)
```

#### RTCUR with resampling Fiber CUR
```
[L_core,X_sub_mat,timer,err ] = RTCUR_rf(D, R, para)
```

## Input Description
1. D : Inputed tensor. 
1. R : Targeted multilinear rank.
1. (optional) para.max_iter, para.epsilon, para.zeta, para.gamma, para.CI: parameters described in our paper. All have defalut values.

* See paper for the details of constant selection.

## Output Description
1. L_core : Core tensor, i.e., $\mathcal{R}$.
1. X_sub_mat : Fiber CUR components, i.e., {$C_i U_i^\dagger$}.

#### To obtain the full estimated tensor, call 
```
L_est = ttm(L_core,X_sub_mat);
```

## Demo

Clone the codes as well as <a href=https://gitlab.com/tensors/tensor_toolbox/-/releases/v3.1>Tensor Toolbox v3.1</a> and run `simpleRPCAtest.m` for a test demo.
