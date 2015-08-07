
import sys

def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print >> sys.stderr, 'fix_pareto.py in.json out.json'
    s = open(args[0], 'rt').read()
    s = s.replace('Infinity', '100')
    with open(args[1], 'wt') as f:
        f.write(s)
        
if __name__ == '__main__':
    main()

