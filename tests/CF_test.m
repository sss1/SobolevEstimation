% This is the original CF test proposed by Epps and Singleton (1986).
function p = CF_test(Xs, Ys, K, J)

  n = size(Xs, 1); % sample size

  D = size(Xs, 2); % data dimension

  % randomly choose IID standard normally test frequencies
  Zs = normrnd(0, 1, D, J);

  Ws = exp(-i * Xs * Zs) - exp(-i * Ys * Zs);
  W = mean(Ws, 1);
  Sigma = cov(Ws);
  T = n*real(W*inv(Sigma)*W'); % Hotelling's T^2 statistic

  df = size(Zs, 2); % degrees of freedom

  p = chi2cdf(T, df, 'upper');

end
