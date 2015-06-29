function [fm, fs, xtrain, ytrain] = gpucb(x, f, fgp, sn, niter, h, plot)
  rng(42);
  xtrain = zeros(0, 1);
  ytrain = zeros(0, 1);
  [~, ~, fm, fs] = fgp(xtrain, ytrain, x);
  noise = randn(niter, 1);
  rng(142);
  for t = 1:niter
    ut = fm + 2*sqrt(fs);
    umax = max(ut);
    imax = find(ut == umax);
    idx = imax(randsample(length(imax), 1));
    xtrain = [xtrain; x(idx)];
    ytrain = [ytrain; f(idx) + sn*noise(t)];
    [~, ~, fm, fs] = fgp(xtrain, ytrain, x);
  end
  if plot
    figure;
    hold on;
    ylim([-2, 3.5]);
    plot(x, f, '--');
    hp = h * ones(length(x), 1);
    plot(x, hp, 'k--');
    plot(x, fm, 'r-');
    plot(x, fm + 2*sqrt(fs), 'r--');
    plot(x, fm - 2*sqrt(fs), 'r--');
    plot(xtrain, ytrain, 'ro');
  end
end