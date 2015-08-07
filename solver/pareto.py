
from collect_pareto import get_pareto
import sys

def main():
    args = sys.argv[1:]
    
    if len(args) == 0:
        print >> sys.stderr, 'python pareto.py in.txt'
        print >> sys.stderr, 'Assumes in.txt has lines with format "time error ...", prints lines on Pareto frontier'
    
    lines = open(args[0], 'rU').read().strip().split('\n')
    lines_split = [x.split(' ', 2) for x in lines]
    time_L = [float(x[0]) for x in lines_split]
    error_L = [float(x[1]) for x in lines_split]
    idx_L = get_pareto(error_L, time_L)

    for i in idx_L:
        print lines[i]
        
if __name__ == '__main__':
    main()
    
