"""

    Define a target function for evaluating the sampler
"""
import numpy
import collections
from scipy.stats import norm


def make_function(N=10, dp=0.1, sd_scale=5, domain=100):
    """
    N - how many gaussians
    dp -- dirichlet parameter. As you go to 0, gets peakier
    domain - [0,domain] is what the function is defined on
    sd_scale -- the scale of the exponential that we sample the sd from
    """

    mu = numpy.random.random(N)*domain  # the mean of each gaussian
    sd = numpy.random.exponential(sd_scale, size=N)
    w  = numpy.random.dirichlet([dp]*N, size=N)  # the weight assigned to each gaussian

    def f_(x):
        assert not isinstance(x, collections.Iterable)  # NOT arrays
        if x < 0 or x > 100:
            return 0
        return numpy.sum(w * norm.pdf((x-mu)/sd))

    def cdf_(x):
        if x < 0:
            return 0
        if x > 100:
            x = 100
        return numpy.sum(w * norm.cdf((x-mu)/sd)) - numpy.sum(w * norm.cdf((0-mu)/sd))

    # Find the true normalizing constant
    normalizingConastant = numpy.sum(w * norm.cdf((100-mu)/sd)) - numpy.sum(w * norm.cdf((0-mu)/sd))

    # Return a function and a cdf
    return lambda x: f_(x) / normalizingConastant, lambda x: cdf_(x) / normalizingConastant


if __name__ == "__main__":

    x = numpy.arange(0, 100, 0.01)
    f = make_function()

    import matplotlib.pyplot as plt

    plt.plot(x, map(f,x), '-')
    plt.show()
