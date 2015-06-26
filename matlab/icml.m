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

%s0 = 100;
%plot(x(s0), f(s0), 'go');

%L = lip(x, f);
%plot(x, f(s0) + L*(-x + x(s0)), 'b');
%plot(x, f(s0) + L*(x - x(s0)), 'b');

%itrain = [30, 130, 160];
%xtrain = x(itrain);
%ytrain = f(itrain) + sn*gpml_randn(0.2, length(itrain), 1);
%[~, ~, fm, fs] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, xtrain, ytrain, x);

%plot(x, fm, 'r-');
%plot(x, fm + 2*sqrt(fs), 'r--');
%plot(x, fm - 2*sqrt(fs), 'r--');
%plot(xtrain, ytrain, 'ro');

%figure;
%niter = 100;
%fgp = @(xtrain, ytrain, x)gp(hyp, @infExact, meanfunc, covfunc, likfunc, xtrain, ytrain, x);
%gpucb(x, f, fgp, sn, niter, h);

figure;
h = 0;
L = lip(x, f);
s0 = 70;
niter = 35;
fgp = @(xtrain, ytrain, x)gp(hyp, @infExact, meanfunc, covfunc, likfunc, xtrain, ytrain, x);
safeopt(x, f, fgp, sn, niter, h, L, s0);