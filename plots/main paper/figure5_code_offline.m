clc
clear all

% Problem setting
n = 2;
m = 1;
sys.A = [0.9802, 0.0196; 0, 0.9802];
sys.C = [1, -1];
Q = [1.9608, 0.0195; 0.0195, 1.9605];
R = 1;
sys.B = [Q^0.5, zeros(2,1)];
sys.D = [0, 0, 1];

run_count = 500;
T = 1000;

% Learn the Kalman gain
x_0 = zeros(2,1);
V_0 = eye(2);
all_rho = (1:0.1:2)*1e-1;
all_c = (1:0.1:2)*1e-4;
G_w = zeros(n,m,T,length(all_rho));
for k = 1 : length(all_rho)
    [~, ~, G_w(:,:,:,k)] = WKF(sys, all_rho(k), nan(m, T), x_0, V_0);
end
[~, ~, G_kalman] = WKF(sys, 0, nan(m, T), x_0, V_0);

%%
% Set seed number
rng(12345);

% Initialization
err_KF = zeros(T, run_count);
err_WKF = zeros(T, run_count, length(all_c));
err_KLKF = zeros(T, run_count, length(all_c));

% Parameters
is_TV = true;   % false
coeff = 0.99;  % 0.99
tau = 0;

for r = 1 : run_count

    fprintf('Runinig Iteration %d\n', r);

    [x, y, y0] = generate_data(sys, x_0, T, coeff, is_TV);

    xhat_kalman = apply_kalman_gain(sys, G_kalman, y, x_0);
    err_KF(:,r) = sum((x-xhat_kalman).^2,1);

    for k = 1 : length(all_rho)
        xhat = apply_kalman_gain(sys, G_w(:,:,:,k), y, x_0);
        err_WKF(:,r,k) = sum((x-xhat).^2,1);
    end

    y_delay = [y0, y(:,1:end-1)]';
    for k = 1 : length(all_c)
        xhat = rkalman(sys, all_c(k), tau, y_delay, x_0, V_0);
        err_KLKF(:,r,k) = sum((x-xhat').^2,1);
    end
end


tmp = mean(mean(err_WKF,2),1); tmp = tmp(:);
[~,k_rho] = min(tmp);
tmp = mean(mean(err_KLKF,2),1); tmp = tmp(:);
[~,k_c] = min(tmp);
figure;
smt = 20;
semilogx(smooth(10*log10(mean(err_KF, 2)), smt), 'LineWidth', 2); hold on;
semilogx(smooth(10*log10(mean(err_KLKF(:,:,k_c), 2)), smt), 'LineWidth', 2);
semilogx(smooth(10*log10(mean(err_WKF(:,:,k_rho), 2)), smt), 'LineWidth', 2);
font_size = 24;
set(gca, 'FontSize', 18);
ylabel('Empirical error (dB)','FontSize', font_size, 'Interpreter', 'latex');
xlabel('time step','FontSize', font_size, 'Interpreter', 'latex');
leg1 = legend({'Kalman filter', 'Kullback-Leibler filter', 'Wasserstein filter'}, 'Interpreter', 'latex', 'Location', 'southeast');
set(leg1,'FontSize',15);
remove_border()

fname = 'large_';
if coeff == 0.099
    fname = 'small_';
end
if is_TV
    fname = strcat(fname, 'tv');
else
    fname = strcat(fname, 'fix');
end

cd figs/
saveas(gcf, fname, 'epsc')
saveas(gcf, fname, 'svg')
cd ..