# Clustering using Symmetric Non-negative Matrix Factorization (syNMF)


## Introduction 

Symmetric nonnegative matrix factorization (SymNMF) is an unsupervised algorithm for graph clustering, 
In NMF we implicity choose our similarity measure as inner products, a fact that can harm the ability to differentiate between different clusters, in contrast, symNMF is a bit more flexible in the ways we choose similiarities for the data points.

## The Algorithm


Given a set of n points 
$X = x_1, x_2, . . . , x_N \in R^d$

the algorithm is:

![Figure1](https://github.com/EladShaba/SymmNMF/blob/main/SymNMF%20algorithm.jpg)

For a more complete understanding of the algorithm, we highly recommend reading the original paper, that goes into much more detalis and possible ways to implement the algorithm [research paper by Da Kuang, Chris Ding and Haesun Park](https://faculty.cc.gatech.edu/~hpark/papers/DaDingParkSDM12.pdf).

## Summary of the functions

API:

Additional files:
* `first.c`:
* `second.c`:

## Basic usage 

## Results 


## Comparing to Kmeans++

Link to github implementation


## Additional information

