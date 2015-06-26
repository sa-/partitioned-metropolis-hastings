from MHSampler import getAbsoluteDifference
from TargetFunction import make_function
import numpy
import sys


## Plotting for varyable distance and steps
for i in range(1000):
    for partitionCount in [10, 100, 1000, 10000]:

        steps = 10000
        partitions = (numpy.random.random(partitionCount) * 100).tolist()
        partitions.append(0)
        partitions.append(100)
        partitions.sort()

        f, cdf   = make_function()

        classic   = getAbsoluteDifference(f, cdf, numberOfIterations=steps, plot=True)
        partition = getAbsoluteDifference(f, cdf, partitions=partitions, numberOfIterations=steps, plot=True)

        print partitionCount, i, classic, "\t", partition
        sys.stdout.flush()
