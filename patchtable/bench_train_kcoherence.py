
import os

def system(s):
  os.system(s)

#os.system('rm bench_train_kcoherence.txt')

for prop_iters in [1, 2, 3, 4, 5]: #[1, 2, 3]:
  for kcoherence in [6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19]: #[1, 2, 3, 4, 5, 10, 20]:
    system('python bench.py ours 1 "-kcoherence %(kcoherence)d -limit 0.01 -init_random 1 -prop_iters %(prop_iters)d -query_step 1 -do_table_lookup 0 -verbose 1 -spatial 0 -rs_max 0 -triangle_factor 1000000000" -one_vs_many 1 >> bench_train_kcoherence.txt' % locals())
