function run_gp(x, f, fgp, sn, itrain)
  if ~isempty(itrain)
    xtrain = x(itrain);
    ytrain = f(itrain) + sn*gpml_randn(0.2, length(itrain), 1);
  else
    xtrain = zeros(0, 1);
    ytrain = zeros(0, 1);
  end
  [~, ~, fm, fs] = fgp(xtrain, ytrain, x);
  prefix = ['gp-', num2str(length(itrain)), '-'];
  csvwrite([prefix, 'mt.csv'], [x, fm]);
  csvwrite([prefix, 'ut.csv'], [x, fm + 2.1*sqrt(fs)]);
  csvwrite([prefix, 'lt.csv'], [x, fm - 2.1*sqrt(fs)]);
  csvwrite([prefix, 'train.csv'], [xtrain ytrain]);
end