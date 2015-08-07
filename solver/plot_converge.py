
import sys
import numpy
import glob
import fnmatch
#import scipy.interpolate
import scipy.spatial
import os.path
from parse_args import parse_args

has_display = 'DISPLAY' in os.environ
if not has_display:
    import matplotlib
    matplotlib.use('Agg')
import pylab
import plot_trace

avg_key = 'avg_time'
max_value = 50.0
#avg_key = 'avg_time20'
#avg_key = 'avg_time10'
#plot_calls = False #True

def avg_time_func(A):
#    print A.shape
    T = A[:,0]
    E = A[:,1]
#    print T.shape
    T = numpy.concatenate(([T[0]], T, [T[-1]]))
    E = numpy.concatenate(([1.0], E, [0.0]))            # x axis
    T = T[::-1]
    E = E[::-1]
    return numpy.trapz(T, E)

omit = '*firh_firv*'

warn_done = False

def trace_xy(L):
    global warn_done
    xvals = numpy.array([(x['of_calls'] if plot_calls else x['wall_time']/60.0) for x in L])
    try:
        yvals = numpy.array([(x[avg_key] if x[avg_key] < 1e6 else max_value) for x in L])
    except KeyError:
        if not warn_done or True:
            print 'Warning: using deprecated integration to find avg_time'
            warn_done = True
        yvals = numpy.array([avg_time_func(x['pareto']) for x in L])
    return (xvals, yvals)
    
def avg_trace(filename0, filenameL):
    traceL = [single_trace(filename) for filename in filenameL]
    XL = [numpy.max(X) for (X, Y) in traceL]
    max_filename = max([(XL[i], filenameL[i]) for i in range(len(filenameL))])[1]
    print 'avg_trace: %-30s (%d), converge: min: %5.2f, mean: %5.2f, median: %5.2f, max: %5.2f (%s)' % (filename0, len(traceL), numpy.min(XL), numpy.mean(XL), numpy.median(XL), numpy.max(XL), max_filename)
    Xmax = 0.0
    for (X, Y) in traceL:
        Xmax = max(Xmax, numpy.max(X))
    Xans = numpy.arange(0.0, Xmax, Xmax/100.0)
    Yans = numpy.zeros(len(Xans))
    for (i, (X, Y)) in enumerate(traceL):
        #print X.shape, X.dtype
        #print Y.shape, Y.dtype
        #print Xans.shape, Xans.dtype
        #print filenameL[i]
        #print 'X'
        #print X
        #print 'Y'
        #print Y
        #print 'Xans'
        #print Xans
        # Workaround for griddata() not properly filling in nearest points (it leaves NaNs), and NearestNDInterpolator not working with 1D data
        #Yinterp = scipy.interpolate.NearestNDInterpolator(numpy.vstack((X, X)).T, Y)(numpy.vstack((Xans, Xans)).T)
        K = scipy.spatial.KDTree(X.reshape((len(X),1)))
        Yinterp = numpy.array([Y[K.query([Xv])[1]] for Xv in Xans])
        #scipy.interpolate.griddata(X, Y, Xans, method='nearest') #, fill_value=0.0)
        #print 'Yinterp'
        #print Yinterp
        #print
        Yans += Yinterp*(1.0/len(filenameL))
    #sys.exit(1)
    return (Xans, Yans)

def single_trace(filename):
    L = plot_trace.read_trace(filename)[0]
    return trace_xy(L)

def usage():
    print >> sys.stderr, 'plot_converge trace1.json [trace2.json ...] [options]'
    print >> sys.stderr, '  If wildcards are used ("*trace.json", quotes needed) then average convergence is plotted'
    print >> sys.stderr, '  (firh_firv are omitted). Can use directory which matches all trace files within.'
    print >> sys.stderr, '   -ground ground.json     -- Plot against ground truth'
    print >> sys.stderr, '   -plot_calls b           -- Plot calls if 1, default 0 (plots wall time).'
    print >> sys.stderr, '   -o plot.png             -- Write plot to file (default if no display available).'
    print >> sys.stderr, '   -indiv b                -- Plot on individual plots (default 0)'
    print >> sys.stderr, '   -xmax x, -ymax y        -- Maximum value for x/y axis (default auto)'
    print >> sys.stderr, '   -legend b'
    sys.exit(1)
    
def main():
    (args, kw) = parse_args('ground plot_calls o indiv legend xmax ymax'.split(), usage)
    if len(args) == 0:
        usage()
    
    has_ground = 'ground' in kw
    if has_ground:
        args = [kw['ground']] + args
    global plot_calls
    plot_calls = bool(int(kw.get('plot_calls', 0)))
    has_out_plot = 'o' in kw
    out_plot = kw.get('o', 'plot.png')
    indiv = int(kw.get('indiv', '0'))
    legend = int(kw.get('legend', '1'))
    xmax_arg = float(kw.get('xmax', '-1'))
    ymax_arg = float(kw.get('ymax', '-1'))
    
    Xmax = 0.0
    for (ifilename, filename) in reversed(list(enumerate(args))):
        filename0 = filename
        #if filename == '-':
        #    continue
        if os.path.isdir(filename):
            filename = os.path.join(filename, '*trace.json')
        if '*' in filename:
            (X, Y) = avg_trace(filename0, [x for x in glob.glob(filename) if not fnmatch.fnmatch(x, omit)])
        else:
            (X, Y) = single_trace(filename)
        #wall_time = [x['wall_time'] for x in L]
        #avg_time = [avg_time_func(x['pareto']) for x in L]
        #(wall_time, avg_time) = trace_xy(L)
        #print X
        #print Y
        Xmax = max(Xmax, numpy.max(X))
        if ifilename == 0 and has_ground:
            Ymin = numpy.min(Y)
            Y = [Ymin, Ymin]
            X = [0.0, Xmax]
            filename = 'Ground Truth'
        pylab.plot(X, Y, '-' if ifilename < 7 else '-+', label=filename)
        pylab.xlabel('OF Calls' if plot_calls else 'Wall Time [min]')
        pylab.ylabel('Mean Program Time')
        if legend:
            pylab.legend()
        if indiv:
            pylab.show()
        #pylab.scatter(wall_time, avg_time, facecolors='none')
    if xmax_arg > 0:
        pylab.xlim(0, xmax_arg)
    if ymax_arg > 0:
        (ymin, ymax) = pylab.ylim()
        pylab.ylim(ymin, ymax_arg)
    if has_out_plot or not has_display:
        pylab.savefig(out_plot, dpi=200)
    else:
        pylab.show()
    
if __name__ == '__main__':
    main()
    