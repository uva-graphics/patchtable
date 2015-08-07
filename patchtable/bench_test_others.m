for n=[10 30 60]
  system('rm -rf ../patchmatch/one_vs_many');
  system(sprintf('cp -r ../patchmatch/one_vs_many%d ../patchmatch/one_vs_many', n));
  fprintf('\n');
  fprintf('TreeCANN %d\n', n);
  bench_treecann_pareto;
  fprintf('\n');
  fprintf('PM %d\n', n);
  bench_pm_pareto;
end
