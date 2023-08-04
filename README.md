# Clustering using Symmetric Non-negative Matrix Factorization (symNMF)


## Introduction 

Symmetric nonnegative matrix factorization (SymNMF) is an unsupervised algorithm for graph clustering, 
In NMF we implicity choose our similarity measure as inner products, a fact that can harm the ability to differentiate between different clusters, in contrast, symNMF is a bit more flexible in the ways we choose similiarities for the data points.

## The Algorithm


Given a set of n points 
$X = x_1, x_2, . . . , x_N \in R^d$

the algorithm is:

![Figure1](https://github.com/EladShaba/SymmNMF/blob/main/SymNMF%20algorithm.jpg)

Calculation of H:
First we randomly initialize H using the values from the interval [0, 2 ∗ sqrt(avg(W) / k)].
We mark this initialized H as $H^{0}$, and we iteratively update $H^{i}$:


$H^{i+1}_{ij} \leftarrow H^{i}_{ij} * (\frac{1}{2} + \frac{1}{2} * \frac{((WH^{i})_{ij})}{(H^{i} (H^{i})^{T} H^{i})_{ij}})$

For a more complete understanding of the algorithm, we highly recommend reading the original paper, that goes into much more detalis and possible ways to implement the algorithm [research paper by Da Kuang, Chris Ding and Haesun Park](https://faculty.cc.gatech.edu/~hpark/papers/DaDingParkSDM12.pdf).

## Summary of the files

Main files:
* `symnmf.py`:
* `symnmf.h`:
* `symnmf.c`:
* `symnmfmodule.c`:
* `setup.py`:
* `Makefile` :
  
Additional files:
* `first.c`:
* `second.c`:

## Basic usage 
add explnation

## Results 


## Comparing to Kmeans++

Link to github implementation


## Additional information

