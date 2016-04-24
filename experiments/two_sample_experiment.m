% A basic experiment demonstrating consistency of the estimator for distinguishing
% two distribution that differ via corruption by additive Gaussian noise, using
% a permutation test. This is rather slow, due to the permutation test.

n_trials = 1000;

f = @(x) 1 + (pi < x); % density function (up to multiplicative constant)
upper_bound = 2; % upper bound on f for monte carlo sampling

ns = round(logspace(1, 6, 30)); % sample sizes to try
sigma = 0.1; % Gaussian noise variance
% n_perms = 1000; % number of permutations in permutation test
alpha = 0.05; % Type I error bound

D = 1; % TODO: generalize this
s = 0; % Sobolev order to estimate

% s2 is the assumed order of Sobolev smoothness. s2 must be greater than s to
% guarantee consistency. s2 = 2s + D/4 is the minimum degree of smoothness
% assumed to guarantee O(1/n) mean squared error
s2 = 2;

Zs = @(n) round(n^(2/(4*s2 + D))); % scaling of Z_n with n

ps = zeros(length(ns), n_trials); % allocate space to save p-values

for n_idx = 1:length(ns)

  n = ns(n_idx)

  is_done = zeros(n_trials, 1);

  tic;

  for trial = 1:n_trials

    Xs = monte_carlo_sample(f, 2, ns(n_idx));
    Ys = monte_carlo_sample(f, 2, ns(n_idx)) + normrnd(0, sigma, ns(n_idx), 1);

    ps(n_idx, trial) = asymptotic_test(Xs, Ys, s, Zs(n));

    % print completion percentage
    is_done(trial) = 1;
    percent = 100*mean(is_done);
    to_print = sprintf('%04.2f%% done after %04.2f seconds', percent, toc);
    for i = 1:(length(to_print) + 4), fprintf('\b'); end
    fprintf('%04.2f%% done after %04.2f seconds', percent, toc);

  end

  [n mean(ps(n_idx, :)) Zs(n)]

end

figure;
hold all;
plot(log10(ns), alpha*ones(size(ns)));
errorbar(log10(ns), mean(ps, 2), std(ps, [], 2));
plot(log10(ns), round(arrayfun(Zs, ns)));
