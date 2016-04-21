% Estimates the s-order Sobolev inner product between the densities from which
% Xs and Ys are IID samples, using the Fourier basis truncated to frequencies
% below Z
%
% N.B.: CURRENTLY, THIS ONLY WORKS FOR D = 1
%
% Inputs:
%   Xs - n-by-D matrix of n IID observations of a D-dimensional variable
%   Ys - n-by-D matrix of n IID observations of a D-dimensional variable
%   s - order of the Sobolev norm to compute
%   Z - maximum frequency at which to truncate the approximation
%
% Outputs:
%   S_hat - estimated sobolev inner product

function S_hat = sobolev_inner_product_est(Xs, Ys, s, Z)

  Zs = -Z:Z;

  % TODO: generalize the outer products Xs*Zs and Ys*Zs to D > 1
  p_hats = mean(exp(-i * Xs * Zs), 1);
  q_hats = mean(exp(-i * Ys * Zs), 1);

  coeffs = abs(Zs).^(2*s);

  S_hat = coeffs*real((p_hats.*conj(q_hats))');

end
