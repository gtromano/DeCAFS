**IMPORTANT**: Older versions of DeCAFS have a major bug that severely affect the computational complexity of the procedure. This was fixed from version 3.3.2. Should you have an older version installed (lower than 3.3.2) please make sure you update your DeCAFS package either through CRAN or GitHub. You can check your version number at the bottom of the documentation page of DeCAFS, via `help("DeCAFS")`. 

**WHAT'S NEW**: In addition to the automatic model selection, we introduced a graphical iterative model selection procedure that aids the user in selecting an appropriate model for a given sequence of observations. This tuning procedure can seriously improve performances under more challenging scenarios. More details can be found by checking the documentation: `help("guidedModelSelection")`. 

---
DeCAFS vignette
---

`DeCAFS` is a `c++` implementation for `R` of the DeCAFS algorithm for performing optimal multiple changepoint detection on detecting the change in mean in presence of autocorrelation or random fluctuations in the data sequence.


## Installation and Requirements

### Installing the package

To install the package from Github: 


```r
# devtools::install_github("gtromano/DeCAFS")
library(DeCAFS)
```


Alternatively one could fork this repository, and: 


```r
# install.packages("DeCAFS", repos = NULL, type = "source")
library(DeCAFS)
```


### Requirements for the installation

The packages requires `Rcpp` with compiler support for the `std` library with the `g++14` standard.


### Bugs and further queries

If any bug should be spotted, or for any information regarding this package, please email the package mantainer: `g` dot `romano` at `lancaster.ac.uk`.

## Introduction


### The model

We model a combination of a radom walk process (also known as standard Brownian motion or Wiener Process) and an AR process. 
Let <img src="/tex/941136d38ca0857891338190d63c3156.svg?invert_in_darkmode&sanitize=true" align=middle width=10.239687149999991pt height=14.611878600000017pt/> be a random vectorm then for <img src="/tex/1029cb1e2fc5675c6163bb23d517888d.svg?invert_in_darkmode&sanitize=true" align=middle width=82.46922914999999pt height=21.18721440000001pt/>, 


<p align="center"><img src="/tex/310e3db7a8c7cbe8cc68ba1c0c16c275.svg?invert_in_darkmode&sanitize=true" align=middle width=83.1867366pt height=12.785402849999999pt/></p>

where

<p align="center"><img src="/tex/6ee3f1cf7dec7817061af80e7ac64276.svg?invert_in_darkmode&sanitize=true" align=middle width=321.73628685pt height=20.50407645pt/></p>
and 
<p align="center"><img src="/tex/a322191f90ce11426bf9910abdb3b0a9.svg?invert_in_darkmode&sanitize=true" align=middle width=220.60486964999998pt height=18.312383099999998pt/></p>

Then, DeCAFS solves the following minimization problem: 

<p align="center"><img src="/tex/149243f90d6bc6ee65f8e24cc0c62b53.svg?invert_in_darkmode&sanitize=true" align=middle width=694.4776096499999pt height=46.74512369999999pt/></p>

Where our <img src="/tex/6afac5d05e4b7176de856343996f9dfe.svg?invert_in_darkmode&sanitize=true" align=middle width=64.52400569999999pt height=26.76175259999998pt/>, <img src="/tex/f5d1cca921c74da95a8d3bc6b49b5b7c.svg?invert_in_darkmode&sanitize=true" align=middle width=64.53039284999998pt height=26.76175259999998pt/> and <img src="/tex/f561bfc183f7551f2335a63fed864e10.svg?invert_in_darkmode&sanitize=true" align=middle width=70.43831354999999pt height=24.65753399999998pt/> is an indicator function..

# Quick Start

This demo shows some of the features present in the `DeCAFS` package. 

Three functions at the moment are present in the package:


|functions          |description                                                             |
|:------------------|:-----------------------------------------------------------------------|
|DeCAFS             |Main function to run the DeCAFS algorithm on a sequence of observations |
|dataRWAR           |Generate a realization of a RW+AR process                               |
|estimateParameters |Estimate the parameters of our model                                    |

At the moment only two functions for data generation and parameter estimation are present, and they all are tailored for the Random Walk. Since l2-FPOP can tackle also other Stochastic Processes, more functions are expected to be added.

### A simple example

We will start generating a Random Walk. The function `dataRWAR` takes in:

- the length of the sequence of observations,
- a poisson parameter regulating the probability of seeing a jump,
- the average magnitude of a change,
- the <img src="/tex/d2207092f6f2646c1ceeb203dfd92d1d.svg?invert_in_darkmode&sanitize=true" align=middle width=16.75048154999999pt height=14.15524440000002pt/> and the <img src="/tex/3f4081ec86e300ae2ce8c2e98ba9a781.svg?invert_in_darkmode&sanitize=true" align=middle width=16.578873299999987pt height=14.15524440000002pt/> parameters,
- the autocorrelation parameter <img src="/tex/f50853d41be7d55874e952eb0d80c53e.svg?invert_in_darkmode&sanitize=true" align=middle width=9.794543549999991pt height=22.831056599999986pt/>.


```r
set.seed(42)
Y = dataRWAR(n = 1e3, poisParam = .01, meanGap = 15, phi = .5, sdEta = 3, sdNu = 1)
y = Y[["y"]]
```

Running DeCAFS is fairly straightforward:


```r
res = DeCAFS(y)
```


We can plot the DeCAFS segmentation (red lines), alongside with our real segmentation (dotted blue lines).

![plot of chunk plot1](figure/plot1-1.png)


## Running the algorithm without estimation
Alternatively, we can also pass all the required parameters in order for it to run.
In this case, since we both have an AR and RW component, we will need to pass down both <img src="/tex/a1e144bab8a124dfabebbec0ea0f6609.svg?invert_in_darkmode&sanitize=true" align=middle width=47.53762529999999pt height=21.18721440000001pt/>, <img src="/tex/99b952f73c1262ce58fe3c0e0992935e.svg?invert_in_darkmode&sanitize=true" align=middle width=47.70921869999999pt height=21.18721440000001pt/> and <img src="/tex/b0c302e6e5edbc86f736ee8872f8e0c8.svg?invert_in_darkmode&sanitize=true" align=middle width=44.49760754999999pt height=22.831056599999986pt/>.


```r
res = DeCAFS(y,  beta = 2 * log(length(y)), modelParam = list(sdEta = 3, sdNu = 1, \phi = .7))
```

```
## Error: <text>:1:84: unexpected input
## 1: res = DeCAFS(y,  beta = 2 * log(length(y)), modelParam = list(sdEta = 3, sdNu = 1, \
##                                                                                        ^
```


### Extreme case: Random Walk

Let's say we now have the <img src="/tex/910282a84e2c5f2f8d376a8ceddbe851.svg?invert_in_darkmode&sanitize=true" align=middle width=53.541747599999994pt height=22.831056599999986pt/>. In this case our model simply becomes a random walk plus noise:

<p align="center"><img src="/tex/310e3db7a8c7cbe8cc68ba1c0c16c275.svg?invert_in_darkmode&sanitize=true" align=middle width=83.1867366pt height=12.785402849999999pt/></p>

Our Algorithm is capable of dealing with this extreme situation:


```r
set.seed(44)
Y = dataRWAR(n = 1e3, poisParam = .01, meanGap = 15, phi = 0, sdEta = 2, sdNu = 1)
y = Y[["y"]]

res = DeCAFS(y,  beta = 2 * log(length(y)), modelParam = list(sdEta = 2, sdNu = 1, phi = 0))
```

which leads to the result:

![plot of chunk plot2](figure/plot2-1.png)


### Extreme case: Autoregressive model

Secondly, let's say that the <img src="/tex/0f5504265f5b5a44d782e0d1fe69fc41.svg?invert_in_darkmode&sanitize=true" align=middle width=47.53762529999999pt height=26.76175259999998pt/> In this case we end up with an Autoregressive model with changes.

In this case we need to set <img src="/tex/643eca8e013cd79579b825eb2e361679.svg?invert_in_darkmode&sanitize=true" align=middle width=47.53762529999999pt height=21.18721440000001pt/>, and for <img src="/tex/2e5f91817369fa1adad8fc24f4787f0f.svg?invert_in_darkmode&sanitize=true" align=middle width=60.93602624999999pt height=22.831056599999986pt/>:


```r
set.seed(46)
Y = dataRWAR(n = 1e3, poisParam = .01, meanGap = 10, phi = .98, sdEta = 0, sdNu = 2)
y = Y[["y"]]

res = DeCAFS(y,  beta = 2 * log(length(y)), modelParam = list(sdEta = 0, sdNu = 2, phi = .98))
```

which leads to the result:

![plot of chunk plot3](figure/plot3-1.png)

we see that in this case we miss one changepoint.

## Contributing to this package

If you have interest to contribute to this package, please do not esitate to contact the maintainer:  `g` dot `romano` at `lancaster.ac.uk`.
