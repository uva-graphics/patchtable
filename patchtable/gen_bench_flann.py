
dim_list = [6, 7, 8, 10, 12]

def main():
    f = open('bench_flann.sh', 'wt')
    arrow = '>'
    for ndims in dim_list:
        for nchroma in [1, 2, 3]:
            for flann_checks in [16, 64, 128]:
                for flann_trees in [4, 8]:
                    f.write('python bench.py ours 4 "-ndims %(ndims)d -nchroma %(nchroma)d -query_iters 1 -flann_checks %(flann_checks)s -flann_trees %(flann_trees)s -lookup_algo kdtree" %(arrow)s bench.txt\n' % locals())
                    arrow = '>>'

if __name__ == '__main__':
    main()
    
                

