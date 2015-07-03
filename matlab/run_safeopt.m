function run_safeopt(x, f, fgp, sn, s0, h, niter, plot)
  L = lip(x, f);
  [fm, fs, xtrain, ytrain, st, mt, gt] = safeopt(x, f, fgp, sn, niter, h, L, s0, plot);
  prefix = ['safeopt-', num2str(niter), '-'];
  csvwrite([prefix, 'ut.csv'], [x, fm + 2.1*sqrt(fs)]);
  csvwrite([prefix, 'lt.csv'], [x, fm - 2.1*sqrt(fs)]);
  csvwrite([prefix, 'train.csv'], [xtrain ytrain]);
  csvwrite([prefix, 'st.csv'], [x(st), -1.88*ones(length(st), 1)]);
  csvwrite([prefix, 'gt.csv'], [x(gt), -1.68*ones(length(gt), 1)]);
  csvwrite([prefix, 'mt.csv'], [x(mt), -1.48*ones(length(mt), 1)]);
end