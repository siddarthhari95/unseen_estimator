Estimation of Unseen Species

CSE 549 Computational Biology Project

implementation of algorithm described in the paper, https://arxiv.org/pdf/1707.03854.pdf

Given a distribution of various biological species over a period of time, this algorithm estimates the number of new species that could be observed extrapolating to an arbitrary time in the future. 

The algorithm does two parts:
1. Estimates the count using Good-Toulmin Estimator method
2. Estimates the multi-population histogram of species by generating the fingerprint data of given set of population.

To run:
1. Install the dependencies by running 'pip install -r requirements.txt'
2. Check the config file and modify any parameters as needed.
3. Run the code using 'python genesis.py'


Weighted Linear Estimator:  U = −∑(i1,...,im) : ∑ij>0(m∏j=1(−tj)ij)φi1,...,imW(i1,...,im).

Generate Histogram implementation: Refer -> https://arxiv.org/pdf/1707.03854.pdf

Estimated Number of species Observed: E[num observed] =∑(over α) H(α)*(1−∏(i 1..m)P) where P = (1−αi)^ni


Sample Output:

Minimized alpha                   : [0.7, 0.5, 0.4, 0.6, 0.7, 0.0, 0.5]

Histogram hCounts:                : [0, 0, 0, 1, 0, 20, 0] 

trueUnseen:                       : 18              

Expected Number                   : 20.0

Weighted Linear Estimator (uCap): : 25

