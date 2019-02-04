import config
import numpy as np
class DistributionGenerator():
    '''
     class to generate various types of distributions
    '''
    count = 0
    def __init__(self):
        self.count = 0

    # write functions for types of distributions

    def uniformDistribution(self, size):
        #pass
        return list(np.random.random_integers(config.lower_bound, config.upper_bound, size))

    def geometricDistribution(self, size):
        return list(np.random.geometric(p=(config.prob_success), size=size))
    
    def dirichletDistribution(self, size):
        param = [1]*size
        prob_dist = np.random.dirichlet(param, 1)
        dirichlet_dist = [int(prob * config.MaxDistributionSize) for prob in prob_dist[0]]
        return dirichlet_dist

    def generateDistributions(self, dist, size):
        '''
         generate synthetic ditributions
        '''
        distribution = []
        if dist == 'uniform':
            distribution = self.uniformDistribution(size)
        elif dist == 'geometric':
            distribution = self.geometricDistribution(size)
        elif dist == 'dirichlet':
            distribution = self.dirichletDistribution(size)
        return distribution

    #can write other distributions if time permits to use them for our experimentation


    