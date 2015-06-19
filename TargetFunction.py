"""

    Define a target function for evaluating the sampler
"""
import numpy
import collections
from scipy.stats import norm
from scipy.integrate import quad as integrate

def make_function(N=10, dp=0.1, sd_scale=5, domain=100):
    """
    N - how many gaussians
    dp -- dirichlet parameter. As you go to 0, gets peakier
    domain - [0,domain] is what the function is defined on
    sd_scale -- the scale of the exponential that we sample the sd from
    """

    mu = numpy.random.random(N)*domain # the mean of each gaussian
    sd = numpy.random.exponential(sd_scale, size=N)
    w  = numpy.random.dirichlet([dp]*N, size=N) # the weight assigned to each gaussian

    def f_(x):
        assert not isinstance(x, collections.Iterable) # NOT arrays
        return numpy.sum(w * norm.pdf((x-mu)/sd))

    # Find the true normalizing constant
    Z, err = integrate(f_, 0, domain)

    # the normalized version
    return lambda x: f_(x) / Z


if __name__ == "__main__":

    x = numpy.arange(0, 100, 0.01)
    f = make_function()

    import matplotlib.pyplot as plt

    plt.plot(x, map(f,x), '-')
    plt.show()


