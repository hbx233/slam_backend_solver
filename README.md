# Implementation of Graph Optimization Solver SLAM Problems

## Introduction

This is a generic framework for graph optimization, especially for solving optimization problem in **S**imultaneous **L**ocalization **A**nd **M**apping. This is very similar to g2o, but I implemented the whole optimization framework from ground up and used some of the modern features in C++11/14 to make the design more succinct. It's not fully optimized and my implementation of sparse matrix solver is not as efficient as some existing packages, it's only just my way to learn something: implement the whole thing from ground up.

## Software Structure

### 1. Graph Optimization Problem 

A graph optimization problem consists of vertics and edges that connect vertices. Every vertex represents a state variable that we want to optimize, every edge represents an observation (noisy) between vertices it connects to. The graph can be thought as a sping connected system, the edges are springs and some of them are stiff (good measurement), some of them are soft (noisy measurement). The goal of graph optimization is find the configuration of vertices that can make the global spring system's enery (error) lowest.

### 2. Linear Solver

The key part of graph optimization problem is solving the normal equation at each step to obtain the update vector, and we need to design an efficient solver to make it possible for large scale optimization problem. 

**Currently I implemented three types of solvers: **

1. Dense linear solver which purely solve the linear system by matrix inversion;
2. Sparse Schur solver which use sparse schur complement
3. Sparse Cholesky solver which use sparse cholesky decomposition

### 3. Nonlinear Least Square Minimizer

Graph optimization also needs a minimizer to control the optimization process and minimize the total cost. Typical nonlinear least square minimizer include Gauss-Newton, Levenberg-Marquardt and Dog-Leg Algorithm. In this work I only implemented Levenberg-Marquardt algorithm.

## Examples

Examples can be found in `slam_backend_solver/app` folder

### 1. Curve Fitting 

The first example is curve fitting. It only contains one vertex (curve parameters) and many edges only connect to one vertex. Since I use variadic template in `BaseEdge` class, so there's no unary edge and binary edge, you only need to specify the set of vertices' types in the template parameter.

The curve for fitting is $y = exp(ax^2 + bx + c)$

Ground truth: `a = 1; b = 2; c = 1`

Parameters before fitting: `a = 0; b = 0; c = 0;`

Parameters after fitting (30 iterations): `a = 1.05215; b = 1.93838; c = 1.0035`

### 2. Sparse Bundle Adjustment



### 3. Pose Graph Optimization 

####Simulated Pose Graph with 100 poses with 386 constrains along a circle

![pgo](./resources/pgo.gif)





