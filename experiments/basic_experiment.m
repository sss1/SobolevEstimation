% A basic experiment demonstrating consistency of the distance estimator for
% distinguishing two samples corrupted by additive Gaussian noise

n_trials = 10;

f = @(x) 1 + (pi < x); % density function (up to multiplicative constant)
upper_bound = 2; % upper bound on f for monte carlo sampling

ns = round(logspace(0, 6, 40)); % sample sizes to try
sigma = 0.1; % Gaussian noise variance

s = 0; % Sobolev order to estimate
s2 = 2; % order of assumed Sobolev smoothness; must be greater than s to converge

Zs = @(n) n^(2/(4*s2 + 1)); % scaling of Z_n with n

S_est = zeros(length(ns), n_trials); % allocate space to save estimates

tic;

for i = 1:length(ns)

  n = ns(i);

  parfor trial = 1:n_trials

    Xs = monte_carlo_sample(f, 2, ns(i));
    Ys = monte_carlo_sample(f, 2, ns(i)) + normrnd(0, sigma, ns(i), 1);

    S_est(i, trial) = real(sobolev_distance_est(Xs, Ys, s, Zs(n)));

  end

end

toc;

figure;
hold all;
plot(log10(ns), zeros(size(ns)));
errorbar(log10(ns), mean(S_est, 2), std(S_est, [], 2));
plot(log10(ns), round(arrayfun(Zs, ns)));
