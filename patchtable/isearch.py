
# -limit 1:
#6 -filter_dims 2,1,5,0,3,10 0.193392
#7 -filter_dims 5,1,3,28,0,6,2 0.188321
#7 -filter_dims 5,1,3,10,0,6,2 0.188808
#8 -filter_dims 2,1,6,3,0,5,19,10 0.188249
#10 -filter_dims 2,6,19,10,3,0,1,15,5,13 0.1869
#11 -filter_dims 13,19,6,1,3,30,10,0,5,28,2 0.187792

import sys
import os
import random
import threadmap
from bench import *

ndims = 7
max_dim = 40
max_iters = 100000

nproc = 8

def f(L):
    L_str = ','.join(str(x) for x in L)
    s = subprocess.check_output('./patchtable match vidpair0/a.png vidpair0/b.png out.pfm -ndims %d -lookup_algo table -dt_algo prop -verbose 1 -run_dt 1 -do_rs 1 -flann_checks 32 -flann_trees 4 -dt_knn 1 -limit 1 -query_iters 1 -filter_dims %s' % (ndims, L_str), shell=True)
    return float_last_line(s, 'mean_dist:')

def main():
    random.seed(0)
    L = [random.randrange(max_dim) for i in range(ndims)] #range(ndims)
    fbest = f(L)

    def get_improvement(i):
        current = list(L)
        not_included = list(set(range(max_dim)) - set(current))
        i = random.randrange(len(current))
        current[i] = random.choice(not_included)
        fprime = f(current)
        if fprime < fbest:
            return (fprime, current)
        return (fbest, L)
    
    for iter in range(max_iters):
        improved_list = threadmap.map(get_improvement, range(nproc), n=nproc, dynamic=True)
        for (fprime, current) in improved_list:
            if fprime < fbest:
                fbest = fprime
                L = current
        print iter, ','.join(str(x) for x in L), fbest 

if __name__ == '__main__':
    main()

