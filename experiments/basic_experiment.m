% A basic experiment demonstrating consistency of the distance estimator for
% distinguishing two samples corrupted by additive Gaussian noise

n_trials = 100;

f = @(x) 1 + (pi < x); % density function (up to multiplicative constant)
upper_bound = 2; % upper bound on f for monte carlo sampling

ns = round(logspace(1, 4, 20)); % sample sizes to try
sigma = 1.0; % Gaussian noise variance

D = 1; % TODO: generalize this
s = 1; % Sobolev order to estimate
% s2 is the assumed order of Sobolev smoothness. s2 must be greater than s to
% guarantee consistency. s2 = 2s + D/4 is the minimum degree of smoothness
% assumed to guarantee O(1/n) mean squared error
s2 = 2;

Zs = @(n) round(n^(2/(4*s2 + 1))); % scaling of Z_n with n

S_est0 = zeros(length(ns), n_trials); % allocate space to save estimates
S_est1 = zeros(length(ns), n_trials); % allocate space to save estimates
S_est_split0 = zeros(length(ns), n_trials); % allocate space to save estimates
S_est_split1 = zeros(length(ns), n_trials); % allocate space to save estimates

for i = 1:length(ns)

  tic;

  n = ns(i);

  parfor trial = 1:n_trials

    Xs = monte_carlo_sample(f, 2, ns(i));
    Ys = monte_carlo_sample(f, 2, ns(i)) + normrnd(0, sigma, ns(i), 1);

    S_est0(i, trial) = real(sobolev_distance_est(Xs, Ys, 0, Zs(n)));
    S_est1(i, trial) = real(sobolev_distance_est(Xs, Ys, 1, Zs(n)));
    S_est_split0(i, trial) = real(sobolev_distance_est_split(Xs, Ys, 0, Zs(n)));
    S_est_split1(i, trial) = real(sobolev_distance_est_split(Xs, Ys, 1, Zs(n)));

  end

  [n toc Zs(n)]

end

% figure;
% hold all;
% plot(log10(ns), zeros(size(ns)));
% errorbar(log10(ns), mean(S_est, 2), std(S_est, [], 2));
% plot(log10(ns), round(arrayfun(Zs, ns)));
