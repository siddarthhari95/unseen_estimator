import random
import config
class Sampler():
    '''
     picks samples from provided distributions
    '''
    #write a function to pick samples from the provided distributions  
    def generateSamples(self, populations, num_samples, sample_size):
        sample_set = []
        for n in range(num_samples):
            size = random.randint(1, sample_size)
            sample = []
            for s in range(size):
                i = random.randint(0, config.num_pops-1)
                j = random.randint(0, config.MaxDistributionSize-1)
                sample.append(populations[i][j])
            sample_set.append(sample)
        return sample_set
