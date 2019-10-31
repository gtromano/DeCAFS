---
title: "l2FPOP vignette"
---



## Installation and Requirements

### Installing the package

To install the package from Github: 


```r
# devtools::install_github("alghul96/l2FPOP")
library(l2FPOP)
```


Alternatively one could fork this repository, and: 


```r
# install.packages("l2FPOP", repos = NULL, type = "source")
library(l2FPOP)
```


### Requirements for the installation

The packages requires `Rcpp` with compiler support for the `std` library with the `g++14` standard.


### Bugs and further queries

If any bug should be spotted, or for any information regarding this package, please email the package mantainer: `g` dot `romano` at `lancaster.ac.uk`.

## Introduction

`l2-fpop` is a `c++` implementation for `R` of the l2-fpop algorithm for performing optimal multiple changepoint detection on some ill-conditioned problems such as detecting a change in mean of the distribution of a Random Walk or on a AR process or both for a stream of univariate data.


### The model

We model a combination of a radom walk process (also known as standard Brownian motion or Wiener Process) and an AR process. 
Let <img src="/tex/941136d38ca0857891338190d63c3156.svg?invert_in_darkmode&sanitize=true" align=middle width=10.239687149999991pt height=14.611878600000017pt/> be a random vectorm then for <img src="/tex/1029cb1e2fc5675c6163bb23d517888d.svg?invert_in_darkmode&sanitize=true" align=middle width=82.46922914999999pt height=21.18721440000001pt/>, 


<p align="center"><img src="/tex/d12fda9822d833a84838c666365b9665.svg?invert_in_darkmode&sanitize=true" align=middle width=185.21077409999998pt height=14.611878599999999pt/></p>

where

<p align="center"><img src="/tex/07ced93a0453aeb28cdb018a6950885f.svg?invert_in_darkmode&sanitize=true" align=middle width=262.37398605pt height=20.50407645pt/></p>
and 
<p align="center"><img src="/tex/9a342b42a7dc1d37586eeb8381326ba4.svg?invert_in_darkmode&sanitize=true" align=middle width=258.04785599999997pt height=18.312383099999998pt/></p>

Then, l2FPOP solves the following minimization problem: 

<p align="center"><img src="/tex/1cef58ad7c14c4c00d8abc205636cf7c.svg?invert_in_darkmode&sanitize=true" align=middle width=834.3748330499999pt height=53.64026084999999pt/></p>

Where our <img src="/tex/6afac5d05e4b7176de856343996f9dfe.svg?invert_in_darkmode&sanitize=true" align=middle width=64.52400569999999pt height=26.76175259999998pt/>, <img src="/tex/f5d1cca921c74da95a8d3bc6b49b5b7c.svg?invert_in_darkmode&sanitize=true" align=middle width=64.53039284999998pt height=26.76175259999998pt/> and <img src="/tex/f561bfc183f7551f2335a63fed864e10.svg?invert_in_darkmode&sanitize=true" align=middle width=70.43831354999999pt height=24.65753399999998pt/> is an indicator function..

# Quick Start

This demo shows some of the features present in the `l2FPOP` package. 

Three functions at the moment are present in the package:


|functions          |description                                                              |
|:------------------|:------------------------------------------------------------------------|
|l2fpop             |Main function to run the l2-FPOP algorithm on a sequence of observations |
|dataRWAR           |Generate a realization of a RW+AR process                                |
|estimateParameters |Estimate the parameters sigma_x, sigma_y of a RW                         |

At the moment only two functions for data generation and parameter estimation are present, and they all are tailored for the Random Walk. Since l2-FPOP can tackle also other Stochastic Processes, more functions are expected to be added.

### The `l2-fpop` function

The `l2-fpop` function takes as input the following arguments:

- `y`: the sequence of observations we want to run the algorithm on;
- `beta`: the penalty for the l0 norm in our minimization;
- `lambda`: the penalty for the first l2 norm;
- `gamma`: the penalty for the second l2 norm;
- `type`: the type of costraint to apply to the recursion. At the moment only the standard change ("std") is implemented.

### A simple example

We will start generating a Random Walk. The function `dataRWAR` takes in:

- the length of the sequence of observations,
- a poisson parameter regulating the probability of seeing a jump,
- the average magnitude of a change,
- the <img src="/tex/d2207092f6f2646c1ceeb203dfd92d1d.svg?invert_in_darkmode&sanitize=true" align=middle width=16.75048154999999pt height=14.15524440000002pt/> and the <img src="/tex/3f4081ec86e300ae2ce8c2e98ba9a781.svg?invert_in_darkmode&sanitize=true" align=middle width=16.578873299999987pt height=14.15524440000002pt/> parameters,
- the autocorrelation parameter <img src="/tex/f50853d41be7d55874e952eb0d80c53e.svg?invert_in_darkmode&sanitize=true" align=middle width=9.794543549999991pt height=22.831056599999986pt/>.


```r
set.seed(42)
Y = dataRWAR(n = 1e3, poisParam = .01, meanGap = 15, phi = .7, sdEta = 2, sdNi = .3)
y = Y[["y"]]
```

Running l2-FPOP is fairly straightforward. We need to pass all the required parameters in order for it to run.

In this case, since we both have an AR and RW component, we will use <img src="/tex/2e236e01a90352dee7e211cb3704d3ee.svg?invert_in_darkmode&sanitize=true" align=middle width=268.7669787pt height=26.76175259999998pt/> and <img src="/tex/b0c302e6e5edbc86f736ee8872f8e0c8.svg?invert_in_darkmode&sanitize=true" align=middle width=44.49760754999999pt height=22.831056599999986pt/>.


```r
res = l2fpop(y,  beta = 2 * log(length(y)), lambda = 1/(2^2), gamma = 1/(.3)^2, phi = 0.7)
```

We plot l2fpop segmentation (red lines), alongside with our real segmentation (dotted blue lines).

![plot of chunk plot1](figure/plot1-1.png)

...aaaand it worked (this time).

## Contributing to this package

If you have interest to contribute to this package, please do not esitate to contact the maintainer:  `g` dot `romano` at `lancaster.ac.uk`.
