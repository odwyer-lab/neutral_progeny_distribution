Novelty, popularity, and emergent neutrality: bias in the choice of baby names and lessons for analyzing cultural data
===============
This repository contains open-source code and example data for a project in collaboration with Anne Kandler. 
Specifically, we include here several R functions associated with the results of the paper (arXiv link TBA). The empirical data for each of these functions is a progeny distribution.

## Project Goals

* identify new metrics to test the hypothesis of neutrality in cultural change data
 
* apply these metrics to empirical data

## Overview of the tools in this repository

During this project we derived an analytical expression for the neutral progeny distribution.  The progeny distribution counts the number of individuals "born" to each species/variant type during a given time span.  An example would be to count the number of baby names registered for each distinct name.

### Specific Contents

In this code we introduce three R functions which can be used to reproduce calculations and plots from the paper, as well as to compute neutral parameters using progeny distributions drawn from other sources.

* the function

```{r eval=FALSE}
cumulasym(b,d,K,thresh)
```

will evaluate numerically our expression for the neutral progeny distribution, in the form of a cumulative distribution. The output is the expected fraction of variants with greater than or equal to K individuals. This requires specifying a birth rate, an effective mortality rate, and a value of K, which refers to number of individuals. Note that the expression actually only depends on the ratio of birth to mortality rate. We have an argument "thresh" which is by default = 1, but will allow us to plot the distribution assuming that all variants with fewer than thresh individuals are dropped from the data.

* the function 

```{r eval=FALSE}
nu_MLE(progenybyclass)
```
computes a maximum likelihood estimator for the neutral speciation rate, specifically the combination of parameters 1-b/d. We assume that the progeny distribution is given in the form of a two column matrix, where the first column is a vector of total progeny, and the second is a vector of numbers of species with the number of progeny in the first column

* the function 

```{r eval=FALSE}
plot_prog_vs_neutral(progenybyclass)
```

plots an empirical progeny distribution alongside a neutral distribution (evaluated using its MLE estimator).

* We also include two example progeny distributions, one generated by a simulated neutral process, the other by a simulated process with novelty disadvantage (explained in more detail in our paper).  These are stored in the form of csv files.  To use these, read the csv as

```{r eval=FALSE}
progenybyclass<-read.csv("example_progeny_dbn_1.csv")
```
or
```{r eval=FALSE}
progenybyclass<-read.csv("example_progeny_dbn_2.csv")
```
and use the resulting data to compute ML estimates or plot data using the functions above.

