function safeucb(x, f, fgp, sn, niter, h, L, s0)
  rng(42);
  st = s0;
  xtrain = zeros(0, 1);
  ytrain = zeros(0, 1);
  [~, ~, fm, fs] = fgp(xtrain, ytrain, x);
  noise = randn(niter, 1);
  for t = 1:niter
    lt = fm - 2*sqrt(fs);
    ut = fm + 2*sqrt(fs);
    for i = 1:length(st)
      zi = st(i);
      cert = find(lt(zi) - L*abs(x - x(zi)) >= h);
      st = union(st, cert);
    end
    ust = ut(st);
    umax = max(ust);
    imax = st(ust == umax);
    idx = imax(randsample(length(imax), 1));
    xtrain = [xtrain; x(idx)];
    ytrain = [ytrain; f(idx) + sn*noise(t)];
    [~, ~, fm, fs] = fgp(xtrain, ytrain, x);
  end
  hold on;
  ylim([-2, 3.5]);
  plot(x, f, '--');
  hp = h * ones(length(x), 1);
  plot(x, hp, 'g--');
  plot(x, fm, 'r-');
  plot(x, fm + 2*sqrt(fs), 'r--');
  plot(x, fm - 2*sqrt(fs), 'r--');
  plot(xtrain, ytrain, 'ro');
  plot(x(st), -2, 'bs');
end