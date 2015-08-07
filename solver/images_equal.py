
import Image
import numpy
import sys

def images_equal(filenameA, filenameB, eps=1e-3):
    A = numpy.asarray(Image.open(filenameA), 'float')/255.0
    B = numpy.asarray(Image.open(filenameB), 'float')/255.0
    d = numpy.mean(numpy.abs(A.flatten()-B.flatten()))
    return d < eps
    
def main():
    args = sys.argv[1:]
    if len(args) != 2:
        print >> sys.stderr, 'images_equal.py filenameA filenameB --- Prints 1 or 0'
    print int(images_equal(args[0], args[1]))
    
if __name__ == '__main__':
    main()
    