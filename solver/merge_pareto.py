
import json
import sys

def usage():
    print >> sys.stderr, 'merge_pareto.py pareto_without_tune.json pareto_with_tune.json out.json'
    sys.exit(1)
    
def main():
    args = sys.argv[1:]
    if len(args) < 3:
        usage()
    L = json.loads(open(args[0], 'rt').read())
    L_missing = json.loads(open(args[1], 'rt').read())
    ans = []
    for i in range(len(L)):
        try:
            sub = L[i]
            if not sub.has_key('tune'):
                sub['tune'] = L_missing[i]['tune']
        except KeyError:
            continue
        ans.append(sub)
    print >> sys.stderr, 'Selected %d points' % len(ans)
    
    with open(args[2], 'wt') as f:
        f.write(json.dumps(L))

if __name__ == '__main__':
    main()
