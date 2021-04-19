# bayesian-removal-model
Hierarchical Bayesian removal model used to correct sampling bias and estimate abundance/biomass from removal samples.

This model script is written in R. JAGS is used to implement the Bayesian sampler and must be downloaded prior to running the model (https://sourceforge.net/projects/mcmc-jags/).
For an example, the script focuses on multi-pass fish samples and aims to correct for a size bias in capture probability.

Here are the main steps to using the model:
- Assign all fish to a size class
- Determine the pass catches for every size class-pass combination for each site
- Identify the number of size classes and number of sites
- Convert pass catch data to appropriate data structure
- Input the JAGS sampler parameters
- Write the model file (removal.model.jags)
- Add all necessary data in to a list and identify model parameters to monitor
- Save model output

Data:
The example data consists of fish catch data from three-pass samples from two different sites. Each fish that is captured is measured and weighed. 

Model structure:
The sampling procedure for multi-pass removal sampling involves removing individuals from a closed population from consecutive passes. The catch of size class j at site i is dependent on the population size (n[i,j]) and the probability of capture (q[i,j]). The fish catch in the first pass (c1[i,j]) then follows a binomial distribution where there are two possible outcomes - caught or not caught - and population size represents the number of independent trials. Catches in the subsequent passes (second and third pass) are based on the intial populiaton and the number of individuals removed in the previous passes. After every pass, the number of individuals removed from the population is subtracted from the initial population to obtain the population size for the subsequent pass (n2[i,j] = n[i,j] - c1[i,j]). We assume that capture probability is constant across passes to simplify the model. 

Biomass is the product of population size and geometric mean weight (bio[i,j] = n[i,j] * gmw[i,j]).

Capture probability for each size class-site combination (q[i,j]) is a logistic function of body size (logit(q[i,j]) <- b0[i] + b.size*size[j]). This function also includes a random site effect (b0[i]) to account for among-site variability in capture probability. 
Note - Site-level effects can also be incoporated into the capture probability function. 

Initial population size is assumed to follow a log-normal distribution (nz[i,j] ~ dlnorm(mu.n[j], tau.n[j])) with hyperparameters (mu.n[y] ~ dunif(0,20); tau.n[y] ~ dunif(0.0000001, 1)). 
The prior for the size effect (b.size) was an uninformative normal distribution (b.size ~ dnorm(0, 0.0001)).
The random site effect is assumed to follow a normal distribution (b0[i] ~ dnorm(b0.mu, b0.tau)) with uninformative priors (b0.mu ~ dnorm(0, 0.0001); b0.tau ~ dunif(0,5)). 



