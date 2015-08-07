
import sys

def tobinary(x, digits):
    return bin(x)[2:].zfill(digits)

def hammingdist(a, b, digits):
    a = tobinary(a, digits)
    b = tobinary(b, digits)
    return sum([a[i] != b[i] for i in range(len(a))])

def hammingdist2d((ax, ay), (bx, by), digits):
    a = ax+ay*(2**digits)
    b = bx+by*(2**digits)
    return hammingdist(a, b, 2*digits)
    
def main():
    args = sys.argv[1:]
    if len(args) < 1:
        print >> sys.stderr, 'python hammingseq.py max_digits [stopLen=8]'
        sys.exit(1)
    max_digits = int(args[0])
    stop = 8
    if len(args) > 1:
        stop = int(args[1])
    
    xL = []
    yL = []
    
    for digits in range(0, max_digits+1):
        ans = [(0, 0)]
        visited = set(ans)
        while len(ans) < 4**digits:
            chosen_dist = -1
            chosen = None
            for i in range(2**digits):
                for j in range(2**digits):
                    if (i, j) in visited:
                        continue
                    dist_L = [hammingdist2d((i, j), x, digits) for x in ans]
                    min_dist = min(dist_L)
                    if min_dist > chosen_dist:
                        chosen_dist = min_dist
                        chosen = (i, j)
            if stop is not None and len(ans) >= stop:
                break
            #print 'chosen:', chosen
            #print 'chosen_dist:', chosen_dist
            #print 'visited:', visited
            #print 'ans:', ans
            #print
            ans.append(chosen)
            visited.add(chosen)
        #print digits, ans
        xL.append([x[0] for x in ans])
        yL.append([x[1] for x in ans])
    #print xL
    #print yL
    print 'vector<vector<int> > antialias_subsample_x(%s);' % repr(xL).replace('[', '{').replace(']', '}')
    print 'vector<vector<int> > antialias_subsample_y(%s);' % repr(yL).replace('[', '{').replace(']', '}')
    
if __name__ == '__main__':
    main()