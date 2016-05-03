%This script computes the L2 norm of certain distributions
%repeated times
n_repeated  = 100;

%sample sizes
total_order_size = 5;
orders = 1:total_order_size;
n_samples = 10.^orders;

%Gaussians derivative L2 norms
gaussian_1 = gaussian_derivative_square(0,1);
true_gaussian_derivative_square = integral(gaussian_1,-inf,inf);
estimated_gaussian_derivative_square = zeros(n_repeated,total_order_size);
s = 1;
s2= 2;
Zs = @(n) n^(2/(4*s2 + 1)); % scaling of Z_n with n

for j = 1:n_repeated
    for i = 1:total_order_size
        Xs = randn(n_samples(i),1);
        estimated_gaussian_derivative_square(j,i) = sobolev_norm_est(Xs, s, Zs(n_samples(i)))/2/pi;
    end
end


%plot
errorbar_gaussian_derivative_square = quantile(estimated_gaussian_derivative_square,3);
true_gaussian_derivative_square_multiples = repmat(true_gaussian_derivative_square,[1,total_order_size]);
figure;
hold on;
errorbar(n_samples,errorbar_gaussian_derivative_square(2,:),errorbar_gaussian_derivative_square(2,:)-errorbar_gaussian_derivative_square(1,:),...
    errorbar_gaussian_derivative_square(3,:)-errorbar_gaussian_derivative_square(2,:),'-kd','LineWidth',2);
plot(n_samples,true_gaussian_derivative_square_multiples,'-mx','LineWidth',2);
h_legends = legend('Estimated Distance','True Distance','Location','northeast');
set(h_legends,'FontSize',20);
xlabel('number of samples','FontSize',20);
ylabel('L_2^2','FontSize',20);
set(gca,'xscale','log');
xlim([n_samples(1),n_samples(end)]);
set(gca,'FontSize',20);

saveas(gcf,'./plots/norm_gaussian_derivative_L2','fig');
export_fig(gcf,'./plots/norm_gaussian_derivative_L2.pdf');