meanfunc = {@meanConst}; hyp.mean = 1.1;
covfunc = {@covSEiso}; ell = 1; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.06; hyp.lik = log(sn);
 
n = 200;
x = linspace(0, 10, n)';
K = feval(covfunc{:}, hyp.cov, x) + 0.00001*eye(n);
mu = feval(meanfunc{:}, hyp.mean, x);
R = chol(K);
seed1 = 0.13;
seed2 = 0.2;
f = R'*gpml_randn(seed1, n, 1) + mu;
y = R'*gpml_randn(seed1, n, 1) + mu + exp(hyp.lik)*gpml_randn(0.2, n, 1);
fgp = @(xtrain, ytrain, x)gp(hyp, @infExact, meanfunc, covfunc, likfunc, xtrain, ytrain, x);

h = 0;
s0 = 70;
hp = h * ones(length(x), 1);
csvwrite('ftrue.csv', [x f]);
xsafe = x(f > h);
xunsafe = x(f <= h);
csvwrite('fsafe.csv', [xsafe, -1.88*ones(length(xsafe), 1)]);
csvwrite('funsafe.csv', [xunsafe, -1.88*ones(length(xunsafe), 1)]);
csvwrite('hp.csv', [x hp]);
csvwrite('s0.csv', [x(s0) f(s0)]);

run_gp(x, f, fgp, sn, []);
run_gp(x, f, fgp, sn, 90);
run_gp(x, f, fgp, sn, [50, 90]);
run_gp(x, f, fgp, sn, [50, 90, 170]);

s0 = 95;
L = lip(x, f);
l1 = f(s0) + L*(-x + x(s0));
l2 = f(s0) + L*(x - x(s0));
csvwrite('cert1-s0.csv', [x(s0) f(s0)]);
csvwrite('cert1-lip1.csv', [x, l1]);
csvwrite('cert1-lip2.csv', [x, l2]);
st = find(l1 >=0 & l2 >= 0);
xst = x(st);
csvwrite('cert1-st.csv', [xst, -1.88*ones(length(xst), 1)]);

for k = 3:5
  prefix = ['cert1-', num2str(k)];
  csvwrite([prefix, '-f.csv'], [x(st) f(st)]);
  newst = st;
  for i = 1:length(st)
    zi = st(i);
    cert = find(f(zi) - L*abs(x - x(zi)) >= h);
    newst = union(newst, cert);
  end
  st = newst;
  xst = x(st);
  csvwrite([prefix, '-st.csv'], [xst, -1.88*ones(length(xst), 1)]);
end

prefix = 'cert1-6';
st = 50:138;
csvwrite([prefix, '-f.csv'], [x(st) f(st)]);
xst = x(st);
csvwrite([prefix, '-st.csv'], [xst, -1.88*ones(length(xst), 1)]);

s0 = 95;
L = lip(x, f);
hp = h * ones(length(x), 1);
ep = 0.15;
l11 = f(s0) + ep + L*(-x + x(s0));
l12 = f(s0) + ep + L*(x - x(s0));
l1 = l11 .* (x <= x(s0)) + l12 .* (x > x(s0));
l21 = f(s0) - ep + L*(-x + x(s0));
l22 = f(s0) - ep + L*(x - x(s0));
l2 = l21 .* (x > x(s0)) + l22 .* (x <= x(s0));
csvwrite('cert2-s0.csv', [[x(s0); x(s0)] [f(s0)-ep; f(s0)+ep]]);
csvwrite('cert2-lip1.csv', [x, l1]);
csvwrite('cert2-lip2.csv', [x, l2]);
st = find(l2 >= 0);
xst = x(st);
csvwrite('cert2-st.csv', [xst, -1.88*ones(length(xst), 1)]);

for k = 3:5
  prefix = ['cert', num2str(k)];
  csvwrite([prefix, '-lt.csv'], [x(st) f(st)-ep]);
  csvwrite([prefix, '-ut.csv'], [x(st) f(st)+ep]);
  newst = st;
  for i = 1:length(st)
    zi = st(i);
    cert = find(f(zi) - ep - L*abs(x - x(zi)) >= h);
    newst = union(newst, cert);
  end
  st = newst;
  xst = x(st);
  csvwrite([prefix, '-st.csv'], [xst, -1.88*ones(length(xst), 1)]);
end

prefix = 'cert6';
st = 54:137;
csvwrite([prefix, '-lt.csv'], [x(st) f(st)-ep]);
csvwrite([prefix, '-ut.csv'], [x(st) f(st)+ep]);
xst = x(st);
csvwrite([prefix, '-st.csv'], [xst, -1.88*ones(length(xst), 1)]);

s0 = 70;

fplot = 0;
run_gpucb(x, f, fgp, sn, h, 0, fplot);
run_gpucb(x, f, fgp, sn, h, 5, fplot);
run_gpucb(x, f, fgp, sn, h, 10, fplot);
run_gpucb(x, f, fgp, sn, h, 20, fplot);

fplot = 1;
run_safeucb(x, f, fgp, sn, s0, h, 0, fplot);
run_safeucb(x, f, fgp, sn, s0, h, 5, fplot);
run_safeucb(x, f, fgp, sn, s0, h, 10, fplot);
run_safeucb(x, f, fgp, sn, s0, h, 20, fplot);
run_safeucb(x, f, fgp, sn, s0, h, 50, fplot);

fplot = 0;
run_safeopt(x, f, fgp, sn, s0, h, 0, fplot);
run_safeopt(x, f, fgp, sn, s0, h, 1, fplot);
run_safeopt(x, f, fgp, sn, s0, h, 5, fplot);
run_safeopt(x, f, fgp, sn, s0, h, 10, fplot);
run_safeopt(x, f, fgp, sn, s0, h, 20, fplot);
run_safeopt(x, f, fgp, sn, s0, h, 30, fplot);
run_safeopt(x, f, fgp, sn, s0, h, 35, fplot);
run_safeopt(x, f, fgp, sn, s0, h, 40, fplot);
run_safeopt(x, f, fgp, sn, s0, h, 50, fplot);
run_safeopt(x, f, fgp, sn, s0, h, 100, fplot);