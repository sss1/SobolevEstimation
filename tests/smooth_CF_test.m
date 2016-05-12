% This is the ``smooth characteristic function'' test proposed by
% Chwialkowski, Ramdas, Sejdinovic, and Gretton (2015).
function p = smooth_CF_test(Xs, Ys, J)

  n = size(Xs, 1); % sample size

  D = size(Xs, 2); % data dimension

  Zs = normrnd(0, 1, D, J);

  coeffs_Xs = mvnpdf(Xs);
  coeffs_Ys = mvnpdf(Ys);

  smoothed_CF_Xs = bsxfun(@times, coeffs_Xs, exp(-i * Xs * Zs));
  smoothed_CF_Ys = bsxfun(@times, coeffs_Ys, exp(-i * Ys * Zs));
  Ws = smoothed_CF_Xs - smoothed_CF_Ys;
  
  W = mean(Ws, 1);
  Sigma = cov(Ws);
  T = n*real(W*inv(Sigma)*W'); % Hotelling's T^2 statistic

  df = size(Zs, 2); % degrees of freedom; also equal to (2*Z)^D

  p = chi2cdf(T, df, 'upper');

end
