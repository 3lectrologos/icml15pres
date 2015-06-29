meanfunc = {@meanConst}; hyp.mean = 1.1;
covfunc = {@covSEiso}; ell = 1; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.05; hyp.lik = log(sn);
 
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
csvwrite('hp.csv', [x hp]);
csvwrite('s0.csv', [x(s0) f(s0)]);

run_gp(x, f, fgp, sn, []);
run_gp(x, f, fgp, sn, 90);
run_gp(x, f, fgp, sn, [50, 90]);
run_gp(x, f, fgp, sn, [50, 90, 170]);

%s0 = 100;
%plot(x(s0), f(s0), 'go');

%L = lip(x, f);
%plot(x, f(s0) + L*(-x + x(s0)), 'b');
%plot(x, f(s0) + L*(x - x(s0)), 'b');

%plot(x, fm, 'r-');
%plot(x, fm + 2*sqrt(fs), 'r--');
%plot(x, fm - 2*sqrt(fs), 'r--');
%plot(xtrain, ytrain, 'ro');

%figure;
%niter = 100;
%fgp = @(xtrain, ytrain, x)gp(hyp, @infExact, meanfunc, covfunc, likfunc, xtrain, ytrain, x);
%gpucb(x, f, fgp, sn, niter, h);

plot = 0;
run_gpucb(x, f, fgp, sn, h, 0, plot);
run_gpucb(x, f, fgp, sn, h, 5, plot);
run_gpucb(x, f, fgp, sn, h, 10, plot);
run_gpucb(x, f, fgp, sn, h, 20, plot);

plot = 0;
run_safeopt(x, f, fgp, sn, s0, h, 0, plot);
run_safeopt(x, f, fgp, sn, s0, h, 1, plot);
run_safeopt(x, f, fgp, sn, s0, h, 5, plot);
run_safeopt(x, f, fgp, sn, s0, h, 10, plot);
run_safeopt(x, f, fgp, sn, s0, h, 20, plot);
run_safeopt(x, f, fgp, sn, s0, h, 30, plot);
run_safeopt(x, f, fgp, sn, s0, h, 35, plot);
run_safeopt(x, f, fgp, sn, s0, h, 40, plot);
run_safeopt(x, f, fgp, sn, s0, h, 50, plot);
run_safeopt(x, f, fgp, sn, s0, h, 100, plot);