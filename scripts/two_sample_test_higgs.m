%load data
load('./data/higgs_200000.mat');
%sample sizes
n_repeated = 1000;
D = 4;
n_samples = [1000,2000,4000,8000,12000];

asympt_rejects = zeros(n_repeated, length(n_samples));
cf_rejects = zeros(n_repeated, length(n_samples));
smooth_cf_rejects = zeros(n_repeated,length(n_samples));
test_power_asympt = zeros(length(n_samples),1);
test_power_cf = zeros(length(n_samples),1);
test_power_smooth_cf = zeros(length(n_samples),1);

n_total = min(size(background_features,1),size(higgs_features,1));
s = 0;
s2= 2;

Zs = @(n) n^(2/(4*s2 + D));
Z  = 2;
J = 2;
alpha = 0.05;

for i = 1:n_repeated
    disp(i);
    parfor j = 1:length(n_samples)
       selected_indices = randperm(n_total,n_samples(j)); 
       background_features_selected = background_features(selected_indices,:);
       higgs_features_selected = higgs_features(selected_indices,:);
       p = asymptotic_test(background_features_selected, ...
           higgs_features_selected, s, Z);
       if p < alpha
           asympt_rejects(i,j) = 1;
       end
       
       p = CF_test(background_features_selected, higgs_features_selected, J);
       if p < alpha
           cf_rejects(i,j) = 1;
       end
       
       p = smooth_CF_test(background_features_selected,higgs_features_selected,J)
       if p < alpha
           smooth_cf_rejects(i,j) = 1;
       end
       
    end
end

for i = 1:length(n_samples)
    test_power_asympt(i,:) = nnz(asympt_rejects(:,i))/n_repeated;
    test_power_cf(i,:) = nnz(cf_rejects(:,i))/n_repeated;
    test_power_smooth_cf(i,:) = nnz(smooth_cf_rejects(:,i))/n_repeated;
end



save('./results/higgs_two_sample.mat');

figure;
hold on;

plot(n_samples,test_power_asympt,'-mx','LineWidth',2);
plot(n_samples,test_power_cf,'-b*','LineWidth',2);
plot(n_samples,test_power_smooth_cf,'-ko','LineWidth',2);
h_legends = legend('H^0 Distance Test','CF','Smooth CF','Location','northwest');
set(h_legends,'FontSize',20);
xlabel('number of samples','FontSize',20);
ylabel('L_2^2','FontSize',20);
set(gca,'xscale','log');
xlim([n_samples(1),n_samples(end)]);
set(gca,'FontSize',20);
saveas(gcf,'./plots/higgs_two_sample','fig');
export_fig(gcf,'./plots/higgs_two_sample.pdf');

