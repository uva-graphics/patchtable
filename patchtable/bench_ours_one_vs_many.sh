python bench.py ours 1 "-speed 0 -partition_step 16" -one_vs_many 1 >> bench_ours_one_vs_many.txt
python bench.py ours 1 "-speed 1 -partition_step 16" -one_vs_many 1 >> bench_ours_one_vs_many.txt
python bench.py ours 1 "-speed 2 -partition_step 16" -one_vs_many 1 >> bench_ours_one_vs_many.txt
python bench.py ours 1 "-speed 3 -partition_step 16" -one_vs_many 1 >> bench_ours_one_vs_many.txt
python bench.py ours 1 "-speed 4 -partition_step 16" -one_vs_many 1 >> bench_ours_one_vs_many.txt
python bench.py ours 1 "-speed 5 -partition_step 16" -one_vs_many 1 >> bench_ours_one_vs_many.txt
python bench.py ours 1 "-speed 6 -partition_step 16" -one_vs_many 1 >> bench_ours_one_vs_many.txt
python bench.py ours 1 "-speed 7 -partition_step 16" -one_vs_many 1 >> bench_ours_one_vs_many.txt
python bench.py ours 1 "-speed 8 -partition_step 16" -one_vs_many 1 >> bench_ours_one_vs_many.txt
python bench.py ours 1 "-speed 9 -partition_step 16" -one_vs_many 1 >> bench_ours_one_vs_many.txt

matlab -r "bench_treecann"
