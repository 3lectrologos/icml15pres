function run_safeucb(x, f, fgp, sn, s0, h, niter, plot)
  L = lip(x, f);
  [fm, fs, xtrain, ytrain, st] = safeucb(x, f, fgp, sn, niter, h, L, s0, plot);
  prefix = ['safeucb-', num2str(niter), '-'];
  csvwrite([prefix, 'ut.csv'], [x, fm + 2.1*sqrt(fs)]);
  csvwrite([prefix, 'lt.csv'], [x, fm - 2.1*sqrt(fs)]);
  csvwrite([prefix, 'train.csv'], [xtrain ytrain]);
  csvwrite([prefix, 'st.csv'], [x(st), -1.88*ones(length(st), 1)]);
end