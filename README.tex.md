---
title: "DeCAFS vignette"
---




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

`DeCAFS` is a `c++` implementation for `R` of the DeCAFS algorithm for performing optimal multiple changepoint detection on detecting the change in mean in presence of autocorrelation or random fluctuations in the data sequence.


### The model

We model a combination of a radom walk process (also known as standard Brownian motion or Wiener Process) and an AR process. 
Let $\bf y$ be a random vectorm then for $t=1,\ldots,n$, 


$$
y_t = \mu_t + \epsilon_t
$$

where

$$
\mu_t = \mu_{t-1} + \eta_t + \delta_t, \quad \eta_t \sim N(0, \sigma_\eta^2), \ \delta_t \ \in R
$$
and 
$$
\epsilon_t = \phi \epsilon_{t-1} + \nu_t \quad \nu_t \sim N(0, \sigma_\nu^2)
$$

Then, DeCAFS solves the following minimization problem: 

$$
\min_{\substack{\mu_{1:n}\\ \delta_{2:n}}}(
       (1-\phi^2) \gamma(y_1 - \mu_1)^2 } + \sum_{t = 2}^n  \left[ \lambda (\mu_{t} - \mu_{t-1} - \delta_{t}) ^ 2 
         + \gamma \Big((y_t - \mu_t) - \phi (y_{t-1} - \mu_{t-1})\Big) ^ 2 + \beta \ \mathcal{I}_{\delta_t \neq 0}\right])
$$

Where our $\lambda = 1/\sigma_\eta^2$, $\gamma = 1/\sigma_\nu^2$ and $\mathcal{I} \in \{0, 1 \}$ is an indicator function..

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
- the $\sigma_\nu$ and the $\sigma_\eta$ parameters,
- the autocorrelation parameter $\phi$.


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
In this case, since we both have an AR and RW component, we will need to pass down both $\sigma_\eta = 3$, $\sigma_\nu = 1$ and $\phi = .7$.


```r
res = DeCAFS(y,  beta = 2 * log(length(y)), modelParam = list(sdEta = 3, sdNu = 1, \phi = .7))
```

```
## Error: <text>:1:84: unexpected input
## 1: res = DeCAFS(y,  beta = 2 * log(length(y)), modelParam = list(sdEta = 3, sdNu = 1, \
##                                                                                        ^
```


### Extreme case: Random Walk

Let's say we now have the $phi = 0$. In this case our model simply becomes a random walk plus noise:

$$
y_t = \mu_t + \epsilon_t
$$

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

Secondly, let's say that the $\sigma_\eta^2 = 0$ In this case we end up with an Autoregressive model with changes.

In this case we need to set $\sigma_\eta = 0$, and for $\phi = 0.98$:


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
