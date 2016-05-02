%repeated times
n_repeated  = 50;
%sample sizes
total_order_size = 5;
orders = 1:total_order_size;
n_samples = 10.^orders;

%One uniform distribution one triangle distribution
uniform_1 = uniform_dis(0,1);
triangle_1 = triangle_dis(0,1);
square_difference_uni_tri= square_loss(uniform_1,triangle_1);
true_distance_uni_tri = integral(square_difference_uni_tri,-inf,inf);
estimated_distances_uni_tri = zeros(n_repeated,total_order_size);
s = 0;
s2= 1;
Zs = @(n) n^(2/(4*s2 + 1)); % scaling of Z_n with n

for j = 1:n_repeated
    disp(j);
    for i = 1:total_order_size
        Xs = rand(n_samples(i),1);
        Ys = slicesample(0.5,n_samples(i),'pdf',triangle_1);
        estimated_distances_uni_tri(j,i) = real(sobolev_distance_est(Xs, Ys, s, Zs(n_samples(i))))/2/pi;
    end
end


%plot
errorbar_uni_tri = quantile(estimated_distances_uni_tri,3);
true_uni_tri_multiples = repmat(true_distance_uni_tri,[1,total_order_size]);
figure;
hold on;
errorbar(n_samples,errorbar_uni_tri(2,:),errorbar_uni_tri(2,:)-errorbar_uni_tri(1,:),...
    errorbar_uni_tri(3,:)-errorbar_uni_tri(2,:),'-kd','LineWidth',2);
plot(n_samples,true_uni_tri_multiples,'-mx','LineWidth',2);
h_legends = legend('Estimated Distance','True Distance','Location','southeast');
set(h_legends,'FontSize',20);
xlabel('number of samples','FontSize',20);
ylabel('L_2^2','FontSize',20);
set(gca,'xscale','log');
xlim([n_samples(1),n_samples(end)]);
set(gca,'FontSize',20);
saveas(gcf,'./plots/distance_unif_tri','fig');
export_fig(gcf,'./plots/distance_unif_tri.pdf');