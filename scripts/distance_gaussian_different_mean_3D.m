%repeated times
n_repeated  = 100;
%sample sizes
D = 3;
total_order_size = 5;
orders = 1:total_order_size;
n_samples = 10.^orders;

%Two Gaussians with different mean
mu1 = [0;0;0];
sigma1 = diag([1,1,1]);
mu2 = [1;1;1];
sigma2 = diag([1,1,1]);

square_difference = multinormal_dis(mu1,sigma1,mu2,sigma2);

disp('computing true distance...')
true_distance_gaussian_mean = integral3(square_difference,-inf,inf,-inf,inf,-inf,inf);
estimated_distances_gaussian_mean = zeros(n_repeated,total_order_size);
s = 0;
s2= 4;
Zs = @(n) n^(2/(4*s2 + D)); % scaling of Z_n with n

start_time = tic;

for j = 1:n_repeated
    disp(j);
    for i = 1:total_order_size
        
        Xs = randn(n_samples(i),D);
        Ys = randn(n_samples(i),D) + 1;
        estimated_distances_gaussian_mean(j,i) = real(sobolev_distance_est(Xs, Ys, s, Zs(n_samples(i))))/(2*pi)^D;
    end
end

end_time = toc(start_time);

save('./distance_gaussian_different_mean_3D.mat');

%plot
errorbar_gaussian_mean = quantile(estimated_distances_gaussian_mean,3);
true_gaussian_mean_multiples = repmat(true_distance_gaussian_mean,[1,total_order_size]);
figure;
hold on;
errorbar(n_samples,errorbar_gaussian_mean(2,:),errorbar_gaussian_mean(2,:)-errorbar_gaussian_mean(1,:),...
    errorbar_gaussian_mean(3,:)-errorbar_gaussian_mean(2,:),'-kd','LineWidth',2);
plot(n_samples,true_gaussian_mean_multiples,'-mx','LineWidth',2);
h_legends = legend('Estimated Distance','True Distance','Location','northeast');
set(h_legends,'FontSize',20);
xlabel('number of samples','FontSize',20);
ylabel('L_2^2','FontSize',20);
set(gca,'xscale','log');
xlim([n_samples(1),n_samples(end)]);
set(gca,'FontSize',20);
saveas(gcf,'./plots/distance_gaussian_different_mean_3D','fig');
export_fig(gcf,'./plots/distance_gaussian_different_mean_3D.pdf');

