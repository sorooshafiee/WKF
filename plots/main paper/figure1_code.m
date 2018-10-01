clc
clear all

% We consider the linear system
% y = x + w
% where w is the random gaussian noise
sigma_xx = 1;
sigma_ww = 0.1;
Sigma = [sigma_xx, sigma_xx; sigma_xx, sigma_xx+sigma_ww];
x_dim = 1;

all_rho = 0.5:0.5:2;
W = cell(length(all_rho), 1);
i = 1;
for rho = all_rho
    [~, Q] = Frank_Wolfe(zeros(2,1), Sigma, rho, x_dim);
    W{i} = Q.Sigma;
    i = i + 1;
end

all_rho_tau1 = 3:3:9;
KL = cell(length(all_rho_tau1), 1);
tau = 0;
i = 1;
for rho_KL = all_rho_tau1
    [~, KL{i}] = tau_update(Sigma, rho_KL, tau, x_dim);
    i = i+1;
end

all_rho_tau1 = 3:3:9;
tau1 = cell(length(all_rho_tau1), 1);
tau = 1;
i = 1;
for rho_tau1 = all_rho_tau1
    [~, tau1{i}] = tau_update(Sigma, rho_tau1, tau, x_dim);
    i = i+1;
end
fprintf('Finish generating data. Start plotting...')

%%
% Plot ellipsoid representing the error for Wasserstein
figure;
hold on; axis equal;
plotErrorEllipse([0,0], Sigma, 0.9, ',''r'', ''LineWidth'', 2')
plotErrorEllipse([0,0], W{1}, 0.9, ',''b'', ''LineWidth'', 2')
plotErrorEllipse([0,0], W{3}, 0.9, ',''k'', ''LineWidth'', 2')


y = linspace(-4,4);
[xx, yy, xy, xw, ww] = stat_cal(Sigma);
plot(xy/yy*y, y, '--r', 'LineWidth', 2);

[xx, yy, xy, xw, ww] = stat_cal(W{1});
plot(xy/yy*y, y, '--b', 'LineWidth', 2);

[xx, yy, xy, xw, ww] = stat_cal(W{3});
plot(xy/yy*y, y, '--k', 'LineWidth', 2);
limx = 5;
limy = 3;
xlim([-limx, limx]);
ylim([-limy, limy]);
grid on
%title('Ellipsoid: 90-percentile, Line: $\hat{x} = S^\star_{xy} S^\star_{yy} y$', 'Interpreter', 'latex', 'FontSize', 20);
legend({'$\rho =0 $', '$\rho = 0.5$', '$\rho = 1.5$'}, 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'NW')
%legend({'$\tau =0 $', '$\tau = 2.85$', '$\tau = 8.85$'}, 'Interpreter', 'latex', 'FontSize', 20)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
xticks(sort([-limx:2:limx, 0]))
set(gcf, 'Position', [100, 100, 800, 800])
remove_border()
cd figs/
saveas(gcf, 'Wass', 'svg')
saveas(gcf, 'Wass', 'epsc')
cd ..

%%
% Plot ellipsoid representing the error for KL
figure;
hold on; axis equal;
plotErrorEllipse([0,0], Sigma, 0.9, ',''r'', ''LineWidth'', 2')
plotErrorEllipse([0,0], KL{1}, 0.9, ',''b'', ''LineWidth'', 2')
plotErrorEllipse([0,0], KL{3}, 0.9, ',''k'', ''LineWidth'', 2')


y = linspace(-4,4);
[xx, yy, xy, xw, ww] = stat_cal(Sigma);
plot(xy/yy*y, y, '--r', 'LineWidth', 2);

limx = 5;
limy = 3;
xlim([-limx, limx]);
ylim([-limy, limy]);
grid on
%title('Ellipsoid: 90-percentile, Line: $\hat{x} = S^\star_{xy} S^\star_{yy} y$', 'Interpreter', 'latex', 'FontSize', 20);
%legend({'$\rho =0 $', '$\rho = 0.5$', '$\rho = 1.5$'}, 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'NW')
legend({'$\rho =0 $', '$\rho = 3$', '$\rho = 9$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'NW');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
xticks(sort([-limx:2:limx, 0]))
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);

set(gcf, 'Position', [100, 100, 800, 800])
remove_border()
cd figs/
saveas(gcf, 'KL', 'svg')
saveas(gcf, 'KL', 'epsc')
cd ..


%%
% Plot ellipsoid representing the error for tau1
figure;
hold on; axis equal;
plotErrorEllipse([0,0], Sigma, 0.9, ',''r'', ''LineWidth'', 2')
plotErrorEllipse([0,0], tau1{1}, 0.9, ',''b'', ''LineWidth'', 2')
plotErrorEllipse([0,0], tau1{3}, 0.9, ',''k'', ''LineWidth'', 2')


y = linspace(-4,4);
[xx, yy, xy, xw, ww] = stat_cal(Sigma);
plot(xy/yy*y, y, '--r', 'LineWidth', 2);

limx = 5;
limy = 3;
xlim([-limx, limx]);
ylim([-limy, limy]);
grid on
%title('Ellipsoid: 90-percentile, Line: $\hat{x} = S^\star_{xy} S^\star_{yy} y$', 'Interpreter', 'latex', 'FontSize', 20);
%legend({'$\rho =0 $', '$\rho = 0.5$', '$\rho = 1.5$'}, 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'NW')
legend({'$\rho =0 $', '$\rho = 3$', '$\rho = 9$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'NW');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
xticks(sort([-limx:2:limx, 0]))
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);

set(gcf, 'Position', [100, 100, 800, 800])
remove_border()
cd figs/
saveas(gcf, 'tau1', 'svg')
saveas(gcf, 'tau1', 'epsc')
cd ..