function run_gpucb(x, f, fgp, sn, h, niter, plot)
  [fm, fs, xtrain, ytrain] = gpucb(x, f, fgp, sn, niter, h, plot);
  prefix = ['gpucb-', num2str(niter), '-'];
  csvwrite([prefix, 'ut.csv'], [x, fm + 2.1*sqrt(fs)]);
  csvwrite([prefix, 'lt.csv'], [x, fm - 2.1*sqrt(fs)]);
  csvwrite([prefix, 'train.csv'], [xtrain ytrain]);
end