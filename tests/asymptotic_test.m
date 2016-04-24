% This test is based on Hotelling's T^2 statistic, which has a theoretical
% asymptotic distribution of chi^2(2*Z), for any fixed Z
function p = asymptotic_test(Xs, Ys, s, Z)

  n = size(Xs, 1);

  Zs = [-Z:-1 1:Z];

  coeffs = repmat(abs(Zs).^s, n, 1);

  Vs = coeffs .* (exp(-i * Xs * Zs) - exp(-i * Ys * Zs));

  W = mean(Vs, 1);
  Sigma = cov(Vs);

  T = n*real(W*inv(Sigma)*W'); % Hotelling's T^2 statistic

  df = length(Zs); % degrees of freedom

  p = chi2cdf(T, df, 'upper');

end
