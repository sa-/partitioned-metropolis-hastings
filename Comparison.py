from MHSampler import getAbsoluteDifference
from TargetFunction import make_function
from scipy.integrate import quad as integrate
import numpy
import sys



## Plotting for varyable distance and steps
for i in range(1000):
    for partitionCount in [10, 100, 1000, 10000]:

        steps = 1000
        partitions = (numpy.random.random(partitionCount) * 100).tolist()
        partitions.sort()

        f   = make_function()
        cdf = lambda x: integrate(f, 0, x)

        classic   = getAbsoluteDifference(f, cdf, numberOfIterations=steps)  # poorly named variables
        partition = getAbsoluteDifference(f, cdf, partitions=partitions, numberOfIterations=steps)

        print partitionCount, i, classic, "\t", partition
        sys.stdout.flush()


# Plot number of partitions
#for npartitions in [1,5,10,100]:
