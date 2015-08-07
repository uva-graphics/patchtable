
import sys
import pylab
import numpy
import random

def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print >> sys.stderr, 'vis_kmeans.py points.txt labels.txt [centers.txt] [knn.txt]'
        sys.exit(1)

    points = numpy.loadtxt(args[0])
    labels = numpy.asarray(numpy.loadtxt(args[1]), 'int')
    dims = points.shape[1]
    assert dims == 2

    has_centers = False
    if len(args) > 2:
        has_centers = True
        centers = numpy.loadtxt(args[2])

    has_knn = False
    if len(args) > 3:
        has_knn = True
        knn = numpy.loadtxt(args[3])

    nlabels = int(numpy.amax(labels))+1

    #labels = numpy.vstack((labels,)*dims).T

    step = 10

    points = points[::step,:]
    labels = labels[::step]
    #labels = labels[::step,:]

    colors = [[random.randrange(256) for j in range(3)] for i in range(nlabels)]
    light_colors = [[255 - (255-c)/2 for c in color_v] for color_v in colors]

    def fmt_color(L):
        return '#' + ''.join('%02x' % x for x in L)

    for i in range(nlabels):
        #print (points[labels==i]).shape
        #print (points[labels==i]).shape
        #print
        pylab.scatter(points[labels==i,0], points[labels==i,1], color=fmt_color(light_colors[i]), alpha=0.5)

    if has_knn:
        for i in range(nlabels):
            for k in range(knn.shape[1]):
                pylab.plot([centers[i,0], centers[knn[i,k],0]],
                           [centers[i,1], centers[knn[i,k],1]], color=fmt_color(colors[i]), alpha=0.5)

    if has_centers:
        for i in range(nlabels):
            pylab.scatter(centers[i:i+1,0], centers[i:i+1,1], marker='D', color=fmt_color([c for c in colors[i]]))

    pylab.show()

if __name__ == '__main__':
    main()