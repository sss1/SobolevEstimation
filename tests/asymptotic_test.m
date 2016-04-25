% This test is based on Hotelling's T^2 statistic, which has a theoretical
% asymptotic distribution of chi^2(2*Z), for any fixed Z
function p = asymptotic_test(Xs, Ys, s, Z)

  n = size(Xs, 1);

  Zs = [-Z:-1 1:Z];

  Ws = bsxfun(@times, abs(Zs).^s, exp(-i * Xs * Zs) - exp(-i * Ys * Zs));

  W = mean(Ws, 1);
  Sigma = cov(Ws);

  T = n*real(W*inv(Sigma)*W'); % Hotelling's T^2 statistic

  df = length(Zs); % degrees of freedom

  p = chi2cdf(T, df, 'upper');

end
