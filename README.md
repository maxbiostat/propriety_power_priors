# On the normalised power prior 

Paper: [[Arxiv]](https://arxiv.org/abs/2004.14912)  [[Statistics in Medicine]](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.9124)

A [power prior](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.6728) is a way of constructing an informative prior distribution based on historical data.
The main idea is to raise the likelihood to a certain scalar, `a0`, usually between 0 and 1, in order to control the amount of information borrowed from the historical data. 

In this paper, [Joseph Ibrahim](https://sph.unc.edu/adv_profile/joseph-g-ibrahim-phd/) and [I](https://lmfcarvalho.org/about/) explore the so-called normalised power prior, where `a0` is allowed to vary according a prior distribution.
We prove a few interesting things about the normalising constant `c(a0)` and devise an approximate method to sample from the joint posterior of the tempering scalar `a0`and the parameters of interest (`theta`).

The package with the routines to run the examples in this repository is [npowerPrioR](https://github.com/maxbiostat/npowerPrioR). 
