% Performs a two-sample permutation test using Sobolev distance.
%
% This is called ``naive'' because it could probably be much more
% computationally efficient. Specifically, I think it might be much faster
% to permute the rows of the exponentiated outer product in
% sobolev_inner_product_est.m, rather than repeatedly calling
% sobolev_distance_est (and hence the somewhat slow exp function) many, many
% times.
%
% N.B.: CURRENTLY, THIS ONLY WORKS FOR D = 1
%
% Inputs:
%   Xs - n-by-D matrix of n IID observations of a D-dimensional variable
%   Ys - n-by-D matrix of n IID observations of a D-dimensional variable
%   s - order of the Sobolev norm to compute
%   Z - maximum frequency at which to truncate approximation
%   n_perms - number of permutations to test
%
% Outputs:
%   p - p-value estimated by permutation test

function p = naive_sobolev_two_sample_permutation_test(Xs, Ys, s, Z, n_perms)

  permuted_S_hats = zeros(n_perms, 1);

  S_hat = sobolev_distance_est(Xs, Ys, s, Z);
  joined_data = [Xs; Ys];

  n_X = size(Xs, 1);
  n_Y = size(Ys, 1);

  parfor perm_idx = 1:n_perms
    permuted_data = joined_data(randperm(n_X + n_Y), :);
    permuted_Xs = permuted_data(1:n_X, :);
    permuted_Ys = permuted_data((n_X + 1):end, :);
    permuted_S_hats(perm_idx) = sobolev_distance_est(permuted_Xs, permuted_Ys, s, Z);
  end

  p = mean(S_hat < permuted_S_hats);

end
