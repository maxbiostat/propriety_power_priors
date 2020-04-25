# On the normalised power prior

A [power prior](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.6728) is a way of constructing an informative prior distribution based on historical data.
The main idea is to raise the likelihood to a certain scalar, `a0`, usually between 0 and 1, in order to control the amount of information borrowed from the historical data. 

In this paper, [Joseph Ibrahim](https://sph.unc.edu/adv_profile/joseph-g-ibrahim-phd/) and [I](https://lmfcarvalho.org/about/) explore the so-called normalised power prior, where `a0` is allowed to vary according a prior distribution.
We prove a few interesting things about the normalising constant `c(a0)` and devise an approximate method to sample from the joint posterior of the tempering scalar `a0`and the parameters of interest (`theta`).

You too can use the methods developed here to do your own (normalised) power prior analysis. 
Let's use the Bernoulli experiment in [Neuenschwander et al. (2009)](https://www.ncbi.nlm.nih.gov/pubmed/19735071) as an example.
The steps are:
1. Write your model in Stan, like [this](https://github.com/maxbiostat/propriety_power_priors/blob/master/code/stan/simple_Bernoulli_prior.stan);
2. Then, use [grid_builder.r](https://github.com/maxbiostat/propriety_power_priors/blob/master/code/grid_builder.r) to estimate `c(a0)` at a number of points, for instance as done [here](https://github.com/maxbiostat/propriety_power_priors/blob/master/code/simple_Bernoulli_estimate_c(a0).r);
3. Now you can modify your Stan program to take a dictionary of approximate values for `c(a0)`, like [so](https://github.com/maxbiostat/propriety_power_priors/blob/master/code/stan/simple_Bernoulli_posterior_normalised_approximate.stan);
4. Finally, run your posterior analysis, as exemplified [here](https://github.com/maxbiostat/propriety_power_priors/blob/master/code/simple_Bernoulli_posterior.r);

Any doubts, shoot me a message at `lmax` dot `fgv` at gmail.

Many thanks to Chris Koenig and [Ben Jones](https://www.plymouth.ac.uk/staff/ben-jones) for testing early version of the code.
