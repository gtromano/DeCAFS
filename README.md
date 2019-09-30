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

The packages requires `Rcpp` with compiler support for the `std` library introduced with `g++11`.


### Bugs and further queries

If any bug should be spotted, or for any information regarding this package, please email the package mantainer: `g` dot `romano` at `lancaster.ac.uk`.

## Introduction

`l2-fpop` is a `c++` implementation for `R` of the l2-fpop algorithm for performing optimal multiple changepoint detection on some ill-conditioned problems such as detecting a change in mean of the distribution of a Random Walk or on a AR Model for a stream of univariate data.

The main recursion is a variation from the 2005 algorithm Functional Pruning Optimal Partitioning. Let <img src="/tex/bda932ed515f3d85dc43e487550e10fc.svg?invert_in_darkmode&sanitize=true" align=middle width=121.79969175pt height=22.465723500000017pt/> a series of observations centered on the data generating process <img src="/tex/bd044317942a49071c489f339c0290f9.svg?invert_in_darkmode&sanitize=true" align=middle width=32.018806049999995pt height=22.465723500000017pt/> with underlying signal <img src="/tex/6d0564274ae1a724f8562409b513eae7.svg?invert_in_darkmode&sanitize=true" align=middle width=26.11694744999999pt height=22.831056599999986pt/> where we have no change if <img src="/tex/f8b5cde6bec8f91a244a06a84f0bfeb9.svg?invert_in_darkmode&sanitize=true" align=middle width=64.93146329999999pt height=22.831056599999986pt/> and a change otherwise, and finally let <img src="/tex/5e1d5620dc7bcd919d9143d70d8d0fd3.svg?invert_in_darkmode&sanitize=true" align=middle width=44.984878949999995pt height=24.65753399999998pt/> be a a general statistical loss function, then, l2-fpop finds the optimal solution to the problem:

<p align="center"><img src="/tex/86498c6ce29420c4642493e1fe0f703c.svg?invert_in_darkmode&sanitize=true" align=middle width=469.57916115pt height=27.170795849999998pt/></p>

This is solved via dynamical programming based on the recursion:

<p align="center"><img src="/tex/c47781832279b380c81e11d42507e526.svg?invert_in_darkmode&sanitize=true" align=middle width=412.38927344999996pt height=39.452455349999994pt/></p>

where

<p align="center"><img src="/tex/437a29663e8036be5bd48dbdd77db37d.svg?invert_in_darkmode&sanitize=true" align=middle width=257.77487505pt height=25.7402211pt/></p>

It will follow an example based on the Random Walk. 

### Random Walk

One of the cases on which our algorithm can be applied is the Random Walk, also known as standard Brownian motion or Wiener Process. Let <img src="/tex/feb8fb6e0d2b94a25deda09a72b9a916.svg?invert_in_darkmode&sanitize=true" align=middle width=14.764759349999988pt height=22.55708729999998pt/> be a random vector build in the following way:

<p align="center"><img src="/tex/9a2095a8338ce98da7bdd6706db1431f.svg?invert_in_darkmode&sanitize=true" align=middle width=90.4983486pt height=16.89938415pt/></p>
with

<p align="center"><img src="/tex/8d7f20031b90f94420403f44088b9bf2.svg?invert_in_darkmode&sanitize=true" align=middle width=225.74001779999998pt height=16.438356pt/></p>

where <img src="/tex/c9be722acb577063285f2fbdc4d651b5.svg?invert_in_darkmode&sanitize=true" align=middle width=97.15936064999998pt height=26.76175259999998pt/> and <img src="/tex/ca0201f819bdb6e2c80e9b0b5bf44277.svg?invert_in_darkmode&sanitize=true" align=middle width=97.84615664999998pt height=26.76175259999998pt/>. Then we find a changepoint if <img src="/tex/dfaeef5bfb0e4ba9c40cf60c3e75ecee.svg?invert_in_darkmode&sanitize=true" align=middle width=69.30752729999999pt height=22.831056599999986pt/>. On this framework, our minimization becomes the following: 

<p align="center"><img src="/tex/e9c16705b817ebbbbb23d0a6ba656203.svg?invert_in_darkmode&sanitize=true" align=middle width=605.25544035pt height=66.55531739999999pt/></p>

Where our <img src="/tex/2b16e86b1f5cdcef98896bae27f880cb.svg?invert_in_darkmode&sanitize=true" align=middle width=55.80001844999999pt height=40.47844019999997pt/>.

# Quick Start

This demo shows some of the features present in the `l2FPOP` package. 

Three functions at the moment are present in the package:


|functions          |description                                                              |
|:------------------|:------------------------------------------------------------------------|
|l2fpop             |Main function to run the l2-FPOP algorithm on a sequence of observations |
|dataRW             |Generate a realization of a RW process plus noise                        |
|estimateParameters |Estimate the parameters sigma_x, sigma_y of a RW                         |

At the moment only two functions for data generation and parameter estimation are present, and they all are tailored for the Random Walk. Since l2-FPOP can tackle also other Stochastic Processes, more functions are expected to be added.

### The `l2-fpop` function

The `l2-fpop` function takes as input the following arguments:

- `y`: the sequence of observations we want to run the algorithm on;
- `l0penalty`: the penalty for the l0 norm in our minimization (what is commonly known as `beta` in the litterature);
- `l2penalty`: the penalty for the l2 norm;
- `type`: the type of costraint to apply to the recursion. At the moment only the standard change ("std") and the isotonic regression ("isotonic") are implemented.

### A simple Random Walk example

We will start generating a Random Walk. The function `dataRW` takes in:

- the length of the sequence of observations,
- a poisson parameter regulating the probability of seeing a jump,
- the average magnitude of a change,
- the <img src="/tex/4ebb7b9f51ef56b286b2b327249caa8c.svg?invert_in_darkmode&sanitize=true" align=middle width=16.84748009999999pt height=14.15524440000002pt/> and the <img src="/tex/5cb4e7e42b9def144355e597c609d550.svg?invert_in_darkmode&sanitize=true" align=middle width=16.472713949999992pt height=14.15524440000002pt/>.


```r
set.seed(42)
Y = dataRW(n = 1e3, poisParam = 0.01, meanGap = 20, sdX = 1, sdY = 1)
y = Y[["y"]]
```

Running l2-FPOP is fairly straightforward. We need to pass the <img src="/tex/ce9b0d1765717c60b7915f2a48951a92.svg?invert_in_darkmode&sanitize=true" align=middle width=16.141629899999987pt height=22.831056599999986pt/> parameter for the penalty (called `l0penalty`) as well as the <img src="/tex/22d952fd172ae91ac1817c8f2b3be088.svg?invert_in_darkmode&sanitize=true" align=middle width=16.141629899999987pt height=22.831056599999986pt/> (called `l2penalty`)
In this case, since it's random walk, we will use <img src="/tex/3d80e8fe26cbc23915cf376a0240f0ba.svg?invert_in_darkmode&sanitize=true" align=middle width=119.69567444999997pt height=26.76175259999998pt/> we discussed in the paper and the <img src="/tex/2b16e86b1f5cdcef98896bae27f880cb.svg?invert_in_darkmode&sanitize=true" align=middle width=55.80001844999999pt height=40.47844019999997pt/>.


```r
res = l2fpop(y, l0penalty = 2 * log(length(y)), l2penalty = 1)
```

We plot our segmentation (red lines), alongside with our real segmentation (dotted blue lines).

![plot of chunk plot1](figure/plot1-1.png)

It seems that in this case the algorithm has missed only one change-point. Let's reduce the average jump size and increase the <img src="/tex/5cb4e7e42b9def144355e597c609d550.svg?invert_in_darkmode&sanitize=true" align=middle width=16.472713949999992pt height=14.15524440000002pt/> noise and repeat the experiment:


```r
set.seed(42)
Y = dataRW(n = 1e3, poisParam = 0.01, meanGap = 10, sdX = 1, sdY = 3)
y = Y[["y"]]

res = l2fpop(y, l0penalty = 2 * (3^2) * log(length(y)), l2penalty = (3^2) / 1)
```

We plot again the results:

![plot of chunk plot2](figure/plot2-1.png)

### An autoregressive example

We consider the autoregressive models with jumps:

<p align="center"><img src="/tex/b19bac3fd78e6e5443cf35825efa768b.svg?invert_in_darkmode&sanitize=true" align=middle width=211.73916389999997pt height=16.438356pt/></p>
with <img src="/tex/df8b6224f6c521e639b1ac40f6f5ef97.svg?invert_in_darkmode&sanitize=true" align=middle width=95.0455473pt height=26.76175259999998pt/>.

We will now generate an AR model. The function `dataAR` was designed to do so. It takes as arguments:

- the length of the sequence of observations,
- a poisson parameter regulating the probability of seeing a jump,
- the average magnitude of a change,
- the starting value <img src="/tex/14adeddbb1889c9aba973ba30e7bce77.svg?invert_in_darkmode&sanitize=true" align=middle width=14.61197759999999pt height=14.15524440000002pt/> for the sequence,
- the <img src="/tex/8cda31ed38c6d59d14ebefa440099572.svg?invert_in_darkmode&sanitize=true" align=middle width=9.98290094999999pt height=14.15524440000002pt/> and the <img src="/tex/11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode&sanitize=true" align=middle width=9.423880949999988pt height=14.15524440000002pt/>.


```r
set.seed(42)
Y = dataAR(n = 1e3, poisParam = 0.01, meanGap = 10, y0 = 10, gamma = .95, sd = 1)
y = Y[["y"]]
```

And again, run the algorithm, specifying the now added gamma parameter:


```r
res = l2fpop(y, l0penalty = 2 * .7^2 * log(length(y)), l2penalty = (.95^2) / (2 * 1^2), gamma = .95)

ggplot(data.frame(t = 1:length(y), y), aes(x = t, y = y)) +
  geom_point() +
  geom_vline(xintercept = res[["changepoints"]], color = 2) +
  geom_vline(xintercept = Y[["cp"]], col = 4,  lty = 3)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

### Isotonic regression

This implementation can also perfom isotonic regression case. Using a different constraint function (<img src="/tex/411b775c9044cccae76c5d3830226c2d.svg?invert_in_darkmode&sanitize=true" align=middle width=46.78169759999999pt height=27.705869399999983pt/>), we are now able to pick only up changes. We change the parameter (type) which defines the type of constraint we're using.


```r
set.seed(43)

Y = dataRW(n = 1e3, poisParam = 0.01, meanGap = 10, sdX = 1, sdY = 1)
y = Y[["y"]]

res = l2fpop(y, 2 * log(length(y)), 1, type = "isotonic")

ggplot(data.frame(t = 1:length(y), y), aes(x = t, y = y)) +
  geom_point() +
  geom_vline(xintercept = res[["changepoints"]], color = 2) +
  geom_vline(xintercept = Y[["cp"]], col = 4,  lty = 3)
```

![plot of chunk isotonic](figure/isotonic-1.png)

## Contributing to this package

If you have interest to contribute to this package, please do not esitate to contact the maintainer:  `g` dot `romano` at `lancaster.ac.uk`.
