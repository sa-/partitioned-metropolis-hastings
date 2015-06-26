from math import *
from random import *
import numpy


class MHChain:
    def __init__(self, n, stepSize, startPosition, lb, rb, f, partitions):
        self.stepSize = stepSize
        self.f = f  # pointer to the function to be sampled from
        self.position = startPosition
        self.leftBound = lb
        self.rightBound = rb
        self.proposal = numpy.random.uniform(-stepSize, stepSize, n)  # Creates proposals.
        self.samples = []
        self.samples.append(self.position)
        self.i = 0  # Chain Iterator
        self.partitions = partitions

        # Create an array to keep track of which partition we have proposed to
        self.proposedToRegion = [0] * (len(partitions) - 1)  # [0,...,0]

    def iterate(self):

        candidate = self.position + self.proposal[self.i]
        self.i = self.i+1

        ratio = self.f(candidate)/self.f(self.position)
        acceptanceProb = min([1., ratio])

        withinBounds = self.leftBound <= candidate and candidate <= self.rightBound

        if uniform(0, 1) < acceptanceProb:
            self.whichRegionWasProposedTo(candidate)
            if(withinBounds):
                self.position = candidate

        self.samples.append(self.position)

    def whichRegionWasProposedTo(self, proposal):
        """
        If [0, x1, x2, ..., 100] are the partitions,
        partition[i],partition[i+1]] are the bounds of a region

        TODO: make this O(log(n))
        """
        # Find which region was propsoed to until the last region
        for i in range(len(self.partitions)-1):
            if(proposal > self.partitions[i] and proposal < self.partitions[i+1]):
                self.proposedToRegion[i] += 1
                return
