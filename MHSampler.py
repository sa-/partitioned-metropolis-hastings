from MHChain import *
import matplotlib.pyplot as plt
from random import uniform

"""
    TODO:
    - 'Absolute Difference' calculations need some cleaning up.
"""


def getAbsoluteDifference(f, cdf, partitions=None, binCount=210, stepSize=1., numberOfIterations=20000, plot=False, verbose=False):
    absoluteDifference = None

    if partitions is None:
        absoluteDifference = classicSampler(f, cdf, binCount, stepSize, numberOfIterations, plot, verbose)
    else:
        absoluteDifference = partitionedSampler(f, cdf, partitions, binCount, stepSize, numberOfIterations, plot, verbose)

    return absoluteDifference


def partitionedSampler(f, cdf, partitions, binCount=210, stepSize=1., numberOfIterations=20000, plot=False, verbose=False):

    """
        TODO: See if making the stepSize a function of the size of the partition makes things better?
    """

    numberOfIterations = numberOfIterations/(len(partitions))

    # Create chains
    chains = []

    # Create sampler chains between partitions
    for i in range(len(partitions)-1):
        lb = partitions[i]
        rb = partitions[i+1]
        pos = uniform(lb, rb)
        chain = MHChain(numberOfIterations, stepSize, pos, lb, rb, f, partitions)
        chains.append(chain)

    # Make them iterate one step at a time.
    for i in range(0, numberOfIterations):
        for chain in chains:
            chain.iterate()

    if verbose:
        print "Counts: "
    # Normalize proposition counts so that they reflect probability values
    for chain in chains:
        if verbose:
            print chain.proposedToRegion
        totalPropositions = 0
        smoothingConstant = 1
        for value in chain.proposedToRegion:
            totalPropositions += value + smoothingConstant
        for i in range(0, len(chain.proposedToRegion)):
            # We're adding 1 to every entry in the matrix before normalizing it
            # It's called smoothing
            chain.proposedToRegion[i] = (float(chain.proposedToRegion[i]) + smoothingConstant)/totalPropositions

    if verbose:
        print "\nStochastic matrix: "
    matrix = []  # Set up the stochastic matrix
    for chain in chains:
        if(verbose):
            print(chain.proposedToRegion)
        matrix.append(chain.proposedToRegion)
    matrix = numpy.matrix(matrix)

    steadyStateStochastic = matrix**1000000
    chainWeights = numpy.array(steadyStateStochastic[0]).flatten()
    if verbose:
        print "\nsteadyStateStochastic:"
        print steadyStateStochastic
        print "\nWeights:"
        print chainWeights

    samples = []
    weights = []
    for i in range(len(chains)):
        for j in range(len(chains[i].samples)):
            samples.append(chains[i].samples[j])
            weights.append(chainWeights[i])

    if plot:
        # plt.subplot(2, 1, 1)
        # for chain in chains:
        #     plt.plot(chain.samples, range(len(chain.samples)))

        # plt.subplot(2, 1, 2)
        plt.grid(True)
        x = range(101)
        y = map(f, x)
        plt.hist(samples, weights=weights, bins=binCount, normed=1)
        plt.plot(x, y, 'r')
        plt.show()

    probabilityMass, edges = numpy.histogram(samples, weights=weights, bins=binCount, normed=1)

    # Compute the CDF of the samples
    sampledCDF = []
    previous = 0
    for i in range(len(probabilityMass)):
        sampledCDF.append(previous + probabilityMass[i])
        previous = sampledCDF[i]
    # Normalize the CDF so it sums to 1
    for i in range(len(sampledCDF)):
        sampledCDF[i] = float(sampledCDF[i])/sampledCDF[-1]

    # Compute the true CDF of curve
    trueCDF = []
    for e in edges[1:]:
        val = cdf(e)
        trueCDF.append(val)

    # The difference in the CDFs at every right edge of the histogram bins
    differenceBetweenCDFs = []
    for i in range(len(trueCDF)):
        differenceBetweenCDFs.append(abs(trueCDF[i]-sampledCDF[i]))

    absoluteDifference = sum(differenceBetweenCDFs)
    return absoluteDifference


def classicSampler(f, cdf, binCount=210, stepSize=5., numberOfIterations=20000, plot=False, verbose=False):
    a = 1.
    x = 0

    samples = []
    samples.append(x)

    proposal = numpy.random.uniform(-a, a, numberOfIterations)

    for i in xrange(1, numberOfIterations):
        candidate = x + proposal[i]
        acceptanceProb = min([1., (f(candidate)/f(x))])
        if numpy.random.uniform(0, 1) < acceptanceProb:
            x = candidate
        samples.append(x)

    if plot:
        plt.subplot(2, 1, 1)
        plt.plot(samples, range(len(samples)))

        plt.subplot(2, 1, 2)
        plt.grid(True)
        x = range(101)
        y = map(f, x)
        plt.hist(samples, bins=binCount, normed=1)
        plt.plot(x, y, 'r')
        plt.show()

    probabilityMass, edges = numpy.histogram(samples, bins=binCount)
    sampledCDF = []
    previous = 0
    for i in range(len(probabilityMass)):
        sampledCDF.append(previous + probabilityMass[i])
        previous = sampledCDF[i]

    # Normalize the CDF so it sums to 1
    for i in range(len(sampledCDF)):
        sampledCDF[i] = float(sampledCDF[i])/sampledCDF[-1]

    # Compute the true CDF of curve
    trueCDF = []
    for e in edges[1:]:
        val = cdf(e)
        trueCDF.append(val)

    # The difference in the CDFs at every right edge of the histogram bins
    differenceBetweenCDFs = []
    for i in range(len(trueCDF)):
        differenceBetweenCDFs.append(abs(trueCDF[i]-sampledCDF[i]))

    absoluteDifference = sum(differenceBetweenCDFs)
    return absoluteDifference
