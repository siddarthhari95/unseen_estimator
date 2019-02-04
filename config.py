#tuning various hyperparameters for testing

#extrapolation factor/s
#there can be different extrapolation factors for different 
# populations, can modify this to handle multiple
T = 1

#max size of distribution
MaxDistributionSize = 1000

#min length of a distribution
MinDistributionSize = 10

#Lower bound for values to be picked for distribution
lower_bound = 1

#Upper bound for values to be picked for distribution
upper_bound = 100

#Probability for success in bernoulli trial of geometric distribution
prob_success = 0.1

#size of sample
max_sample_size = 50

# distribution to generate populations from
#dist = 'uniform'

#number of populations
num_pops = 7

# number of samples
num_samples = 7

#distribution types
# dist_types = ["uniform", "dirichlet", "geometric"]
dist_types = ["uniform"]