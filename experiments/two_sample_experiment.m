% A basic experiment demonstrating consistency of the estimator for distinguishing
% two distribution that differ via corruption by additive Gaussian noise, using
% a permutation test. This is rather slow, due to the permutation test.

n_trials = 20;

f = @(x) 1 + (pi < x); % density function (up to multiplicative constant)
upper_bound = 2; % upper bound on f for monte carlo sampling

ns = round(logspace(1, 5, 30)); % sample sizes to try
sigma = 0.1; % Gaussian noise variance
n_perms = 1000; % number of permutations in permutation test
alpha = 0.05; % Type I error bound

s = 0; % Sobolev order to estimate
s2 = 2; % order of assumed Sobolev smoothness; must be greater than s to converge

Zs = @(n) round(n^(2/(4*s2 + 1))); % scaling of Z_n with n

ps = zeros(length(ns), n_trials); % allocate space to save p-values

tic;

for i = 1:length(ns)

  n = ns(i);

  parfor trial = 1:n_trials

    Xs = monte_carlo_sample(f, 2, ns(i));
    Ys = monte_carlo_sample(f, 2, ns(i)) + normrnd(0, sigma, ns(i), 1);

    ps(i, trial) = naive_sobolev_two_sample_permutation_test(Xs, Ys, s, Zs(n), n_perms);

  end

  [n mean(ps(i, :)) Zs(n)]

end

toc;

figure;
hold all;
plot(log10(ns), alpha*ones(size(ns)));
errorbar(log10(ns), mean(ps, 2), std(ps, [], 2));
plot(log10(ns), round(arrayfun(Zs, ns)));
