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
%   alpha - level of confidence interval (only used if CI is requested)
%
% Outputs:
%   S_hat - estimated sobolev inner product
%   CI - confidence interval based on asymptotic normality of the estimator

function [S_hat, CI] = sobolev_inner_product_est(Xs, Ys, s, Z, alpha)

  D = size(Xs, 2);

  Zs = permn(-Z:Z, D)';

  coeffs = prod(abs(Zs).^s, 1);

  if nargout <= 1 % just give point estimate

    % TODO: generalize the outer products Xs*Zs and Ys*Zs to D > 1
    p_hats = mean(exp(-i * Xs * Zs), 1);
    q_hats = mean(exp(-i * Ys * Zs), 1);

    S_hat = (coeffs.^2)*real((p_hats.*conj(q_hats))');

  else if nargout == 2 % estimate confidence interval as well

    if nargin < 5
      alpha = 0.05; % default confidence level
    end

    Ws = bsxfun(@times, coeffs, exp(-i * Xs * Zs));
    Vs = bsxfun(@times, coeffs, exp(-i * Ys * Zs));

    W = mean(Ws, 1);
    V = mean(Vs, 1);
    Sigma_p = cov(Ws);
    Sigma_q = cov(Vs);

    % note that we could have different sample sizes from each density
    sigma2_p = V*Sigma_p*V'/size(Xs, 1);
    sigma2_q = W*Sigma_q*W'/size(Ys, 1);

    S_hat = dot(W, V);
    sigma_hat = sqrt(sigma2_p + sigma2_q); % pool variances

    CI = norminv([alpha/2 (1-(alpha/2))], S_hat, sigma_hat);

  else

    error('Incorrect number of outputs requested; should be 0, 1, or 2.');

  end

end
