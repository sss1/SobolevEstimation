% Estimates the s-order Sobolev distance between the densities from which Xs
% and Ys are IID samples, using the Fourier basis truncated to frequencies
% below Z
%
% N.B.: CURRENTLY, THIS ONLY WORKS FOR D = 1
%
% Inputs:
%   Xs - n-by-D matrix of n observations of a D-dimensional variable
%   Ys - n-by-D matrix of n observations of a D-dimensional variable
%   s - order of the Sobolev norm to compute
%   Z - maximum frequency at which to truncate the approximation
%
% Outputs:
%   S_hat - estimated sobolev distance

function S_hat = sobolev_distance_est(Xs, Ys, s, Z)

  D = size(Xs, 2);

  Zs = permn(-Z:Z, D)';

  % TODO: generalize the outer products Xs*Zs and Ys*Zs to D > 1
  diff_hats = (mean(exp(-i * Xs * Zs), 1) - mean(exp(-i * Ys * Zs), 1))';

  coeffs = prod(abs(Zs).^(2*s), 1);

  S_hat = coeffs*(real(diff_hats).^2 + imag(diff_hats).^2);

end
