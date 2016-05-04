% This test is based on Hotelling's T^2 statistic, which has a theoretical
% asymptotic distribution of chi^2(2*Z), for any fixed Z
function p = asymptotic_test(Xs, Ys, s, Z)

  n = size(Xs, 1); % sample size

  D = size(Xs, 2); % data dimension

  Zs = permn([-Z:-1 1:Z], D)'; % all test (multidimensional) frequencies

  coeffs = prod(abs(Zs).^s, 1);

  Ws = bsxfun(@times, coeffs, exp(-i * Xs * Zs) - exp(-i * Ys * Zs));
  W = mean(Ws, 1);
  Sigma = cov(Ws);
  T = n*real(W*inv(Sigma)*W'); % Hotelling's T^2 statistic

  df = size(Zs, 2); % degrees of freedom; also equal to (2*Z)^D

  p = chi2cdf(T, df, 'upper');

end
