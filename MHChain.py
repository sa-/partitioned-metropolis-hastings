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
        self.proposedToRegion = [0] * (len(partitions) + 1)  # [0,...,0]

    def iterate(self):

        candidate = self.position + self.proposal[self.i]
        self.i = self.i+1

        acceptanceProb = min([1.,
            (self.f(candidate)/self.f(self.position))])

        withinBounds = self.leftBound <= candidate and candidate <= self.rightBound

        if uniform(0, 1) < acceptanceProb:
            self.whichRegionWasProposedTo(candidate)
            if(withinBounds):
                self.position = candidate

        self.samples.append(self.position)

    def whichRegionWasProposedTo(self, proposal):
        # Region 0 is (-inf,partition[0]].
        # Region i is (partition[i+1], partition[i+2]]
        # Region n+1 is (partition[-1],inf)

        # Find which region was propsoed to until the last region
        for i in range(len(self.partitions))[1:]:
            if(proposal < self.partitions[i]):
                self.proposedToRegion[i] += 1
                return

        # If it is in the last region
        if(proposal >= self.partitions[-1]):
            self.proposedToRegion[len(self.partitions)] += 1
