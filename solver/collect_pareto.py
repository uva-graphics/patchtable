
import sys
import json
import glob
import os.path
import time
import numpy
import collections
from parse_args import parse_args

def get_pareto(error_L, time_L):
    """Get indices of Pareto frontier."""
    ERROR = 0
    TIME = 1
    INDEX = 2
    L = sorted([(error_L[idx], time_L[idx], idx) for idx in range(len(error_L))])
    ans = []
    for i in range(len(L)):
        if len(ans) == 0:
            ans.append(L[i][INDEX])
        else:
            error_prev = error_L[ans[-1]]
            time_prev = time_L[ans[-1]]
            error_current = L[i][ERROR]
            time_current = L[i][TIME]
            assert error_current >= error_prev
            if time_current < time_prev:
                ans.append(L[i][INDEX])
    return ans

def usage():
    print >> sys.stderr, 'python collect_pareto.py dir_from_ga [pareto_full_out.txt] [-error_thresh E] [-pareto_layers n] [-rank_top filename.txt] [-subsample n] [-filespec "spec"]'
    sys.exit(1)

def main():
    (args, kw) = parse_args('error_thresh pareto_layers rank_top subsample filespec'.split(), usage)
    error_thresh = float(kw.get('error_thresh', '0.5'))
    pareto_layers = int(kw.get('pareto_layers', '4'))
    rank_top = kw.get('rank_top', '')
    subsample = int(kw.get('subsample', '1'))
    filespec = kw.get('filespec', '*')
    
    if len(args) < 1:
        usage()
    
    outdir = args[0]
    pareto_full = ''
    if len(args) > 1:
        pareto_full = args[1]
    
    T0 = time.time()
    L = []
    avg_errorL = []
    suffix = 'pareto_full.txt'
    filenameL = glob.glob(os.path.join(outdir, filespec + suffix))
    filename_successL = []
    for filename in filenameL:
        try:
            obj_current = json.loads(open(filename, 'rt').read())
            trace_filename = filename[:len(filename)-len(suffix)] + 'trace.json'
#            print trace_filename
            trace_current = json.loads(open(trace_filename, 'rt').read())
            for obj in obj_current:
                obj['filename'] = filename
            L.append(obj_current)
            avg_errorL.append(trace_current[-1]['avg_time20'])
            filename_successL.append(filename)
        except:
            print 'Warning: could not load %s' % filename
    print 'Loaded %d Pareto frontiers in %f secs' % (len(L), time.time()-T0)
    
    object_L = []
    
    error_L = []
    time_L = []
    index_L = []
    for i in range(len(L)):
        for obj in L[i]:
            error_L.append(obj['error'])
            time_L.append(obj['time'])
            filename = os.path.split(filenameL[i])[1]
            prefix = filename.split('_')[0]
            obj['filename'] = os.path.abspath(filenameL[i])
            obj['name'] = prefix + '_' + obj['name']
            object_L.append(obj)
            index_L.append(len(object_L)-1)
    error_L = numpy.array(error_L)
    time_L = numpy.array(time_L)
    index_L = numpy.array(index_L)
    
    print 'Number of scattered points: %d' % len(error_L)
    
    ans = []
    
    for it in range(pareto_layers):
        paretoL = get_pareto(error_L, time_L)
        for idx in paretoL:
            ans.append(object_L[index_L[idx]])
        error_L = numpy.delete(error_L, paretoL)
        time_L = numpy.delete(time_L, paretoL)
        index_L = numpy.delete(index_L, paretoL)
    
    ans = [obj for obj in ans if obj['error'] <= error_thresh]
    if subsample > 1:
        ans = ans[::subsample]
        
    if len(rank_top):
        collect_count = collections.defaultdict(lambda: 0)
        for obj in ans:
            collect_count[obj['filename']] += 1
        count_filename = sorted([(count, filename) for (filename, count) in collect_count.items()], reverse=True)
        
        with open(rank_top, 'wt') as f:
            for (count, filename) in count_filename:
                print >>f, filename
                

    print 'Number of Pareto points written: %d' % len(ans)
    
    if len(pareto_full):
        with open(pareto_full, 'wt') as f:
            f.write(json.dumps(ans)) #, sort_keys=True, indent=4)) #, separators=(',', ': ')))
    
if __name__ == '__main__':
    main()

