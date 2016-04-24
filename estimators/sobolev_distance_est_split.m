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

  Xs_norm = sobolev_norm_est(Xs, s, Z);
  Ys_norm = sobolev_norm_est(Ys, s, Z);
  inner_product = sobolev_inner_product_est(Xs, Ys, s, Z);

  S_hat = Xs_norm + Ys_norm - 2*inner_product;

end
