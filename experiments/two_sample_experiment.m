% A basic experiment demonstrating consistency of the estimator for distinguishing
% two distribution that differ via corruption by additive Gaussian noise, using
% a permutation test. This is rather slow, due to the permutation test.

n_trials = 100;

f = @(x) 1 + (pi < x); % density function (up to multiplicative constant)
upper_bound = 2; % upper bound on f for monte carlo sampling

ns = round(logspace(3, 5, 20)); % sample sizes to try
sigma = 1; % Gaussian noise variance
% n_perms = 1000; % number of permutations in permutation test
alpha = 0.05; % Type I error bound

D = 2; % data dimension
s = 1; % Sobolev order to estimate

% s2 is the assumed order of Sobolev smoothness. s2 must be greater than s to
% guarantee consistency. s2 = 2s + D/4 is the minimum degree of smoothness
% assumed to guarantee O(1/n) mean squared error
s2 = 3;
Zs = @(n) round(sqrt(log10(n))); % @(n) round(n^(2/(4*s2 + D))); % scaling of Z_n with n
J = 10;

ps = zeros(length(ns), n_trials); % allocate space to save p-values

for n_idx = 1:length(ns)

  n = ns(n_idx);

  is_done = zeros(n_trials, 1);

  tic;

  for trial = 1:n_trials

%     Xs = monte_carlo_sample(f, 2, ns(n_idx));
%     Ys = monte_carlo_sample(f, 2, ns(n_idx)) + normrnd(0, sigma, ns(n_idx), 1);
    Xs = normrnd(0, sigma, ns(n_idx), D);
    Ys = normrnd(0.0, sigma, ns(n_idx), D);

    % ps(n_idx, trial) = asymptotic_test(Xs, Ys, s, Zs(n));
    ps(n_idx, trial) = CF_test(Xs, Ys, J);
    % ps(n_idx, trial) = smooth_CF_test(Xs, Ys, J);

    % print completion percentage
    is_done(trial) = 1;
    percent = 100*mean(is_done);
    to_print = sprintf('%04.2f%% done after %04.2f seconds', percent, toc);
    for i = 1:(length(to_print) + 4), fprintf('\b'); end
    fprintf('%04.2f%% done after %04.2f seconds', percent, toc);

  end

  disp(' ');
  disp(['n: ' num2str(n)]);
  disp(['Z: ' num2str((2*Zs(n))^D)]);
  disp(['Power: ' num2str(mean(ps(n_idx, :) < alpha))]);

end

figure(1);
hold all;
plot(log10(ns), alpha*ones(size(ns)));
errorbar(log10(ns), mean(ps, 2), std(ps, [], 2)./sqrt(n_trials));
% plot(log10(ns), round(arrayfun(Zs, ns)));
rejections = ps < alpha;
errorbar(log10(ns), mean(rejections, 2), std(rejections, [], 2)./sqrt(n_trials));
