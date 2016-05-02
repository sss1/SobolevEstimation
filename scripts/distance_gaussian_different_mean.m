%repeated times
n_repeated  = 100;
%sample sizes
total_order_size = 5;
orders = 1:total_order_size;
n_samples = 10.^orders;

%Two Gaussians with different mean
gaussian_1 = gaussian_dis(0,1);
gaussian_2 = gaussian_dis(1,1);
square_difference = square_loss(gaussian_1,gaussian_2);
true_distance_gaussian_mean = integral(square_difference,-inf,inf);
estimated_distances_gaussian_mean = zeros(n_repeated,total_order_size);
s = 0;
s2= 2;
Zs = @(n) n^(2/(4*s2 + 1)); % scaling of Z_n with n

for j = 1:n_repeated
    for i = 1:total_order_size
        Xs = randn(n_samples(i),1);
        Ys = randn(n_samples(i),1) + 1;
        estimated_distances_gaussian_mean(j,i) = real(sobolev_distance_est(Xs, Ys, s, Zs(n_samples(i))))/2/pi;
    end
end

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
saveas(gcf,'./plots/distance_gaussian_different_mean','fig');
export_fig(gcf,'./plots/distance_gaussian_different_mean.pdf');

