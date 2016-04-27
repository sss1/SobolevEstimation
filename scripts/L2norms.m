%This script computes the L2 norm of certain distributions
%sample sizes
total_order_size = 5;
orders = 1:total_order_size;
n_samples = 10.^orders;

%Gaussians L2 norms
gaussian_1 = gaussian_dis_square(0,1);
true_distance_gaussian_mean = integral(gaussian_1,-inf,inf);
estimated_gaussian_norm_square = zeros(1,total_order_size);
s = 0;
s2= 2;
Zs = @(n) n^(2/(4*s2 + 1)); % scaling of Z_n with n

for i = 1:total_order_size
    Xs = randn(n_samples(i),1);
    estimated_gaussian_norm_square(i) = sobolev_norm_est(Xs, s, Zs(n_samples(i)))/2/pi;
end
disp(estimated_gaussian_norm_square - true_distance_gaussian_mean);