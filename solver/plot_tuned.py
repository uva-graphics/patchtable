import sys
import numpy
import pylab
import json
from plot_pareto import colors

def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print >> sys.stderr, 'plot_tuned.py pareto1.json pareto2.json [...]'
        print >> sys.stderr, '  Plots tuned times sorted by first Pareto tuned times'
        sys.exit(1)

    for (i, filename) in enumerate(args):
        pareto = json.loads(open(filename, 'rt').read())
        pareto = sorted(pareto, key=lambda x: x['name'])
        try:
            tuneT = numpy.array([x['tune']['time'] for x in pareto])
        except KeyError:
            tuneT = numpy.array([float(x['baseline']) for x in pareto])
        if i == 0:
            idx = numpy.argsort(tuneT)
        pylab.scatter(numpy.arange(len(tuneT)), tuneT[idx], label=filename, color=colors[i%len(colors)])
        print filename, numpy.mean(tuneT)
    pylab.legend()
    pylab.show()

if __name__ == '__main__':
    main()

