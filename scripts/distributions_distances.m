%This script compare the true L2 distancebetween two functions and
%estimated L2 distance

%sample sizes
total_order_size = 5;
orders = 1:total_order_size;
n_samples = 10.^orders;

%Two Gaussians with different mean
gaussian_1 = gaussian_dis(0,1);
gaussian_2 = gaussian_dis(1,1);
square_difference = square_loss(gaussian_1,gaussian_2);
true_distance_gaussian_mean = integral(square_difference,-inf,inf);
estimated_distances_gaussian_mean = zeros(1,total_order_size);
s = 0;
s2= 2;
Zs = @(n) n^(2/(4*s2 + 1)); % scaling of Z_n with n

for i = 1:total_order_size
%     Xs = slicesample(0,n_samples(i),'pdf',gaussian_1);
%     Ys = slicesample(0,n_samples(i),'pdf',gaussian_2);
    Xs = randn(n_samples(i),1);
    Ys = randn(n_samples(i),1) + 1;
    estimated_distances_gaussian_mean(i) = real(sobolev_distance_est(Xs, Ys, s, Zs(n_samples(i))))/2/pi;
end
disp(estimated_distances_gaussian_mean - true_distance_gaussian_mean);

%Two Gaussians with different vars
gaussian_1 = gaussian_dis(0,1);
gaussian_2 = gaussian_dis(0,2);
square_difference_gaussian_var = square_loss(gaussian_1,gaussian_2);
true_distance_gaussian_var = integral(square_difference_gaussian_var,-inf,inf);
estimated_distances_gaussian_var = zeros(1,total_order_size);
s = 0;
s2= 2;
Zs = @(n) n^(2/(4*s2 + 1)); % scaling of Z_n with n

for i = 1:total_order_size
%     Xs = slicesample(0,n_samples(i),'pdf',gaussian_1);
%     Ys = slicesample(0,n_samples(i),'pdf',gaussian_2);
    Xs = randn(n_samples(i),1);
    Ys = 2*randn(n_samples(i),1);
    estimated_distances_gaussian_var(i) = real(sobolev_distance_est(Xs, Ys, s, Zs(n_samples(i))))/2/pi;
end
disp(estimated_distances_gaussian_var - true_distance_gaussian_var);

%Two uniform distributions with different range
uniform_1 = uniform_dis(0,1);
uniform_2 = uniform_dis(0.5,1.5);
square_difference_uniform = square_loss(uniform_1,uniform_2);
true_distance_uniform = integral(square_difference_uniform,-inf,inf);
estimated_distances_uniform = zeros(1,total_order_size);
s = 0;
s2= 1;
Zs = @(n) n^(2/(4*s2 + 1)); % scaling of Z_n with n

for i = 1:total_order_size
%     Xs = slicesample(0,n_samples(i),'pdf',gaussian_1);
%     Ys = slicesample(0,n_samples(i),'pdf',gaussian_2);
    Xs = rand(n_samples(i),1);
    Ys = 0.5+rand(n_samples(i),1);
    estimated_distances_uniform(i) = real(sobolev_distance_est(Xs, Ys, s, Zs(n_samples(i))))/2/pi;
end
disp(estimated_distances_uniform - true_distance_uniform);


%One uniform distribution one triangle distribution
uniform_1 = uniform_dis(0,1);
triangle_1 = triangle_dis(0,1);
square_difference_uni_tri= square_loss(uniform_1,triangle_1);
true_distance_uni_tri = integral(square_difference_uni_tri,-inf,inf);
estimated_distances_uni_tri = zeros(1,total_order_size);
s = 0;
s2= 1;
Zs = @(n) n^(2/(4*s2 + 1)); % scaling of Z_n with n

for i = 1:total_order_size
    Xs = rand(n_samples(i),1);
%     Ys = 0.5+rand(n_samples(i),1);
    Ys = slicesample(0.5,n_samples(i),'pdf',triangle_1);
    estimated_distances_uni_tri(i) = real(sobolev_distance_est(Xs, Ys, s, Zs(n_samples(i))))/2/pi;
end
disp(estimated_distances_uni_tri - true_distance_uni_tri);

