
README 
--------
This project contains code for using the weighted simplicial complex unmixing m$

There are several key folders to the project.  The doc folder contains README i$

The src folder contains the core source code. The post processing folder contai$
code to post-process inferences that are not simplicial complexes.  The extras $


Dependencies
-------------
SeDuMi, MVES, 2-stage medoidshift clustering, maximal Cliques, and fuzzy membership clustering.

SeDuMi is an optimization package used by MVES
MVES is minimum volume enclosing simplex code used in previous versions of the work as well as this work.
2-stage medoidshift embeds the ideas presented in Roman et al. 2016 for use here.
getFuzzyMembershipClustering allows for mixed membership in each of the phylogenetic branches of tumor evolution.
maximalCliques finds maximal cliques in an adjacency matrix.

USAGE
--------

[V,A,data,VProj,M,I,E,K,CIp,Ps,Fs,RW,penDf,penPf,dimList] = ...
 weightedSCUnmixMainSliver(data,theta,gamma,c,sigma)
is the main function for the core unmixing problem

Inputs:
-------
data is a samples * features matrix in principal components space of the desired gene / protein features
theta is the parameter used in the sliver-based dimensionality detection module (it is a scalar between 0 and 1).  The larger this value, the higher the perceived noise in the data.
gamma is the parameter used to control the weight of the conditional versus prior portion of the probability model used as the objective function.  Higher values of gamma represent higher noise levels.  Gamma = 0 corresponds to the 'hard unmix' version of the problem presented in Schwartz and Shackney, 2010.  We gamma =1.1 in the corresponding manuscript.
c is the coefficient matrix (loading matrix) relating scores (data) to the original values in gene space (non-PC space).
sigma is the parameter to control for the number of standard deviations past the mean to look for in the cluster-wise dimensionality detection module.  Higher values of sigma correspond to more conservative cluster-wise dimension estimates.  We use sigma =  3 in the corresponding manuscript.

Outputs:
--------
V is a matrix representing inferred vertices correspondent to subpopulation genome profiles (in reduced dimensionality space)
A is the adjacency matrix of the simplicial body where vertices are outlined by V.  It is also a matrix.
data is an echo of the data that was unmixed
VProj is a matrix projecting V into the full-dimensional PC space of data.
M represents the percentage-wise membership of each data point in each cluster/sub-simplex/phylogenetic branch.
I,E,K,CIp,Ps are the eigenvalues, errors, initial vertices, and projections, as outlined in earlier works.  They are presented, but not used in downstream analysis for the corresponding publication.
Fs are the mixture fractions in each sub-simplex for each data point.
RWs are the re-weights -- the recomputed weights for each data point for each vertex.
penPf and penDf are the penalty functions for the prior (Pf) and distance / conditional (Df) portions of the final iteration of the objective function.
dimList is the list of dimensions used for each sub-simplex.

If the returned A does not represent 1 connected component, the post-processing functions must be used.  The main post-processing function is as follows:
[Vnew,Anew]=determine_min_sc(V,A,s)

Input
----- 
V and A are the VProj and A generated from the main core function
s is the representation of data used as input for the main unmixing function

Output
-------
Vnew and Anew are the replacements for VProj and A from the core unmixing.

If the core unmixing is used, the M and F values will need to be recomputed.  The named functions perform those tasks.

This code is copyright 2017 Carnegie Mellon University, and can be freely used and modified.  The code was written by Theodore Roman, under the supervision of advisor Russell Schwartz.  


Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
