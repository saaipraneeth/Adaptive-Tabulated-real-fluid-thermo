# quickhull2d.py (modified)
#
# This is a pure Python version of the Quick Hull algorithm.
# It's based on the version of literateprograms, but fixes some
# old-style Numeric function calls.
#
# This version works with numpy version 1.2.1
#
# See: http://en.literateprograms.org/Quickhull_(Python,_arrays)
#
# Modified by Wim Bakker 20081210
#

import numpy
import matplotlib.pyplot as plt
def qhull(sample):
    link = lambda a,b: numpy.concatenate((a,b[1:]))
    edge = lambda a,b: numpy.concatenate(([a],[b]))

    def dome(sample,base): 
        h, t = base
        dists = numpy.dot(sample-h, numpy.dot(((0,-1),(1,0)),(t-h)))
        outer = numpy.repeat(sample, dists>0, axis=0)
        
        if len(outer):
            pivot = sample[numpy.argmax(dists)]
            return link(dome(outer, edge(h, pivot)),
                        dome(outer, edge(pivot, t)))
        else:
            return base

    if len(sample) > 2:
        axis = sample[:,0]
        base = numpy.take(sample, [numpy.argmin(axis), numpy.argmax(axis)], axis=0)
        return link(dome(sample, base),
                    dome(sample, base[::-1]))
    else:
        return sample



# MAIN
if __name__ == "__main__":
    from pylab import plot

    sample = 100*numpy.random.random((32,2))
    hull = qhull(sample)
    plt.figure() 
    for s in sample:
        plt.plot([s[0]], [s[1]], 'b.')

    i = 0
    #while i < len(hull)-1:
    #    plt.plot([hull[i][0], hull[i+1][0]], [hull[i][1], hull[i+1][1]], color='k')
    #    i = i + 1

    #plt.plot([hull[-1][0], hull[0][0]], [hull[-1][1], hull[0][1]], color='k')
    plt.savefig('test.png')
