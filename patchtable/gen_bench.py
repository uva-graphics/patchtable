
import random
#min_dim = 8
#max_dim = 12   # Was 8
#dim_step = 2
#dim_list = [6, 7, 8, 10, 12]
dim_list = [6, 7, 8, 9, 10] #[10, 6, 7, 8, 12]

"""
    for ndims in range(min_dim, max_dim+1):
        for limit in [0.1, 0.5, 1.0, 4.0]:
            for (do_rs, do_prop) in [(0, 0), (0, 1), (1, 1)]:
                for run_dt in [0, 1]:
                    for partition_step in [16, 32]:
                        for dt_iters in [-1, 1]:
                            for query_step in [1, 2, 3]:
"""

def main():
    random.seed(0)
    L = []
    for ndims in dim_list:
        for nchroma in [1]: #, 3]:
            for limit in [1.0, 10.0, 100.0]:
                for kcoherence in [0, 1, 2, 3, 4, 5, 10, 20]:                       # Was [..., 20]
                    for prop_iters in [1, 2, 3]:                                    # Was [1, 2]
                        for spatial in [0, 1]:                                      # Was [0, 1, 2]
                            for (do_rs, do_prop) in [(0, 1)]:                       # Was [(0, 1), (1, 1)]
                                for run_dt in [1]: #0, 1]:
                                    for partition_step in [1]:
                                        for dt_iters in [-1]: #[-1, 1]:
                                            for query_step in [1, 2, 3, 4, 5, 6, 10]: #, 2]:
                                                for kcoherence_iter in [-1]:        # Was [-1, 0]:
                                                    for ntables in range(1, 7):     # Was range(1, 7):
                                                        for triangle_factor in [1.0, 1.5, 2.0, 4.0]:    # Was not searched previously
#                                            for treecann in [0, 1]:
#                                                treecann_paramL = [(0,0)] if not treecann else [(x, y) for x in range(1,11,2) for y in range(1,11,2)]
#                                                for (treecann_agrid, treecann_bgrid) in treecann_paramL:
#-treecann %(treecann)d -treecann_agrid %(treecann_agrid)d -treecann_bgrid %(treecann_bgrid)d
                                                            L.append('python bench.py ours 1 "-ndims %(ndims)d -nchroma %(nchroma)d -limit %(limit)f -do_rs %(do_rs)d -threads 1 -run_dt %(run_dt)d -partition_step %(partition_step)d -dt_iters %(dt_iters)d -do_prop %(do_prop)d -query_step %(query_step)d -kcoherence %(kcoherence)d -kcoherence_iter %(kcoherence_iter)d -prop_iters %(prop_iters)d -spatial %(spatial)d -ntables %(ntables)d -triangle_factor %(triangle_factor)f" -one_vs_many 1 >> bench.txt' % locals())
    
    random.shuffle(L)
    f = open('bench.sh', 'wt')
    L[0] = L[0].replace('>>', '>')
    f.write('\n'.join(L))

if __name__ == '__main__':
    main()

