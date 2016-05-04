% Estimates the s-order Sobolev norm of the density from which Xs is an IID
% sample, using the Fourier basis truncated to frequencies below Z
%
% Inputs:
%   Xs - n-by-D matrix of n observations of a D-dimensional variable
%   s - order of the Sobolev norm to compute
%   Z - maximum frequency at which to truncate the approximation
%
% Outputs:
%   S_hat - estimated sobolev norm

function [S_hat, CI] = sobolev_norm_est(Xs, s, Z)

  n = floor(size(Xs, 1)/2);

  Xs_1 = Xs(1:n, :);
  Xs_2 = Xs((n + 1):end, :);

  [S_hat, CI] = sobolev_inner_product_est(Xs_1, Xs_2, s, Z);

  S_hat = real(S_hat);
  CI = real(CI);

end
