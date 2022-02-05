## Feb. 4
### [Pilanci, 2017](https://stanford.edu/~pilanci/papers/PilWai17.pdf)

* Paper proposes Newton sketch which generalizes IHS to general constrained and unconstrained convex optimization problems. One difference is that in least squares problems, the unsketched Hessian is the same across iterations.
* For d > n problems, a dual strategy exists that provides same guarantees as n > d problems.
* For linear programming problems with n constraints and d dimensions, standard Newton's method for barrier method is O(nd^2) whereas with ROS sketching is O(ndlog(m)+md^2). Sketching is cheaper for n >> d. Partial sketching strategies for additive optimization functions where the Hessians of some subfunctions are sketched and others are retained exactly are also possible.
* Size m lower bound is dependent on a specified tolerance and the Gaussian width of a cone space for the given problem. Euclidean error decays at a linear-quadratic convergence rate.

## Feb. 2
### [Chowdhury, 2018](https://proceedings.mlr.press/v80/chowdhury18a/chowdhury18a.pdf)

* Paper proposes an iterative sketch-based method (with fixed sketches) for ridge regression in the case n << d. A structural constraint in terms of the effective degrees of freedom of the design matrix is imposed that removes the need for the sketch size to be proportional to n. The theory bounds the relative error of the approximate solution to the optimal solution. Appendix results suggest that refreshed matrices help with convergence.

### [Lacotte, 2019](https://proceedings.neurips.cc/paper/2019/file/51425b752a0b402ed3effc83fc4bbb74-Paper.pdf)

* Paper proposes an adaptive sampling method for convex optimization problems where the adaptive sketching matrix also depends on the data matrix. The relative error of the approximate solution w.r.t. the optimal solution is bounded above with high probability and depends on the spectral decay of the data matrix (in addition to other properties). An iterative method is also introduced where the sketching matrix is fixed across iterations.

## Feb. 1
### Yang, 2017

* Paper examines a sketch-based approximation of kernel ridge regression, which is hindered by cubic and squared complexity in terms of number of observations for time and space, respectively. The proposed approach is based on solving the original KRR optimization problem where the empirical kernel matrix is replaced by a sketched version. Their results suggest that if the dimension of the sketched matrix is greater than some value proportional to the statistical dimension of the kernel matrix (the number of eigenvalues greater than some specified tolerance), than the approximate KRR solution still satisfies the minimax optimality (error in sketched function to true function w.r.t. RKHS norm) up to tolerance. These results hold for K-satisfiable matrices, which include Gaussian and randomized orthogonal system sketches (discrete Fourier transform, Hadamard) if the dimension is large enough. TODO: relation to project? Because bounds on error given in terms of sketched solution and true solution and not sketched and original solution?

## Jan. 28

* Skimmed through multiple references (see related literature). There appears to be many works that are very similar but under different titles (presumably previous versions of the same work). Some seem to have been extended while stuck in limbo (arXiv) which leads to confusion about chronological ordering of works.

## Jan. 27
### Dobriban, 2019

* First to be able to distinguish performance of different sketching matrices.
* Considers various loss functions.
* Note: possible connections to factor rotations? See uniform sampling section.

## Jan. 26
### [Pilanci, 2016](https://arxiv.org/pdf/1411.0347.pdf)

* IHS introduced with discussion of least squares problems with convex constraints.
* Original motivation is to address sub-optimality (w.r.t. solution approximation) of classical least squares sketching.
* Sketch dimension per iteration only needs to be proportional to dimension of optimal solution.
* log(1/epsilon) iterations for epsilon-accurate solution.
* Ideas: regularized regression, convex Lp-norm regression problems.

## Jan. 23
### Couillet, 2011

* More comprehensive introduction to random matrix theory.

## Jan. 19-22
### Lacotte, 2020

* The convergence rate results are derived based on the limiting spectral distributions, which are obtained by applying results from random matrix theory/free probability.

### Xia, 2019

* Light introduction and decent starting point to random matrix theory.
* Note: possible connections to factor rotations? See cumulants section.