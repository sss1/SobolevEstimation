%repeated times
n_repeated  = 100;
%sample sizes
total_order_size = 5;
orders = 1:total_order_size;
n_samples = 10.^orders;

%Two uniform distributions with different range
offset = 0.5;
uniform_1 = uniform_dis(0,1);
uniform_2 = uniform_dis(offset,offset+1);
square_difference_uniform = square_loss(uniform_1,uniform_2);
true_distance_uniform_range = integral(square_difference_uniform,-inf,inf);
estimated_distances_uniform = zeros(n_repeated,total_order_size);
s = 0;
s2= 1;
Zs = @(n) n^(2/(4*s2 + 1)); % scaling of Z_n with n

for j = 1:n_repeated
    for i = 1:total_order_size
        Xs = rand(n_samples(i),1);
        Ys = offset+rand(n_samples(i),1);
        estimated_distances_uniform(j,i) = real(sobolev_distance_est(Xs, Ys, s, Zs(n_samples(i))))/2/pi;
    end
end

%plot
errorbar_uniform_range = quantile(estimated_distances_uniform,3);
true_uniform_range_multiples = repmat(true_distance_uniform_range,[1,total_order_size]);
figure;
hold on;
errorbar(n_samples,errorbar_uniform_range(2,:),errorbar_uniform_range(2,:)-errorbar_uniform_range(1,:),...
    errorbar_uniform_range(3,:)-errorbar_uniform_range(2,:),'-kd','LineWidth',2);
plot(n_samples,true_uniform_range_multiples,'-mx','LineWidth',2);
h_legends = legend('Estimated Distance','True Distance','Location','southeast');
set(h_legends,'FontSize',20);
xlabel('number of samples','FontSize',20);
ylabel('L_2^2','FontSize',20);
set(gca,'xscale','log');
xlim([n_samples(1),n_samples(end)]);
set(gca,'FontSize',20);
saveas(gcf,'./plots/distance_unif_different_range','fig');
export_fig(gcf,'./plots/distance_unif_different_range.pdf');