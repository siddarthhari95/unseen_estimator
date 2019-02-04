import random
from generator import DistributionGenerator
import config
import math as mt
import pandas as pd
import scipy.stats as scstat
def generatePopulations(num_pops):
    """
    Generate num_pops populations of the type dist_type.
    """
    distGenerator = DistributionGenerator()
    populations = []
    for i in range(num_pops):
        dist_type = random.sample(config.dist_types, 1)[0]
        populations.append(distGenerator.generateDistributions(dist_type, config.MaxDistributionSize))
    return populations

def generateAlphas(num_pops):
    # list_of_lists = []
    alphas_list = []
    for j in range(num_pops):
        alphas_list.append(round(random.random(),1))
    # list_of_lists.append(alphas_list)
    return alphas_list

def generateHAlphas(alphas, prop_dists):
    # hAlpha = []
    # for j in range(len(alphas)):
    ha = []
    for i in range(len(alphas)):
        if alphas[i] in prop_dists[i]:
            ha.append(prop_dists[i].count(alphas[i]))
        else:
            ha.append(0)
    # hAlpha.append(ha)
    return ha


def generateProbDistributions(pops):
    pop_dists = []
    for pop in pops:
        s = pd.Series(pop)
        pop_dists2 = list((s.groupby(s).transform('count') / len(s)).values)
        pop_dists.append([round(float(i),1) for i in pop_dists2])
    return pop_dists

def poisson(lamba,n):
    # return (mt.exp(lamba) * (lamba**n)) / mt.factorial(n)
    p = mt.exp(-lamba)
    for i in xrange(n):
        p *= lamba
        p /= i+1
    if p == 0:
        p += 1e-6
    return p

def binomial(p, n, i):
    # print p
    p = round(p, 2)
    if p == 1.0:
        p -= 1e-8
    if p == 0.0:
        p += 1e-8
    log = mt.log
    exp = mt.exp
    output = n*log(1-p)
    for i in range(0,i):
        output += log(n-i) + log(p/(1-p)) - log(i+1)
    return exp(output)

def poissonCoefficient(lmd, n):
    log = mt.log
    exp = mt.exp
    output = sum((log(lmd) - log(i+1)) for i in range(0,n))
    return exp(output -lmd)


# if __name__ == "__main__":
#     print(generate_alphas(5))