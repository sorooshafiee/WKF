clc
clearvars

% Set seed number
rng(12345);

% Problem setting
opts.verbose = false;
opts.tol = 1e-4;
run_count = 100;
all_d = 10:10:100;

% Initialization
iter_count = zeros(2, length(all_d), run_count);
time_count = zeros(2, length(all_d), run_count);

for r = 1 : run_count

    fprintf('Runinig Iteration %d\n', r);

    k = 1;
    iter_r = zeros(2, length(all_d));
    time_r = zeros(2, length(all_d));

    for d = all_d

        rho = sqrt(d);
        n = fix(d/5);

        A = randn(d);
        [R, ~] = eig(A + A');
        lambda = 0.1 + 9.9*rand(d,1);
        Sigma = R * diag(lambda) * R';

        t3 = cputime();
        [~, ~, obj] = Frank_Wolfe(zeros(d,1), Sigma, rho, n, opts);
        t4 = cputime();

        time_r(2,k) = t4 - t3;
        iter_r(2,k) = length(obj);
        k = k + 1;

    end

    time_count(:,:,r) = time_r;
    iter_count(:,:,r) = iter_r;
end

%%
figure
shaded = true;
font_size = 18;
h1 = plot(all_d, mean(time_count(2,:,:),3)', 'linewidth', 4);
grid on
set(gca, 'FontSize', 12);
ylabel('Execution time (s)','FontSize', font_size, 'Interpreter', 'latex');
xlabel('Dimension $(d)$','FontSize', font_size, 'Interpreter', 'latex');
if shaded
    prc = 0;
    alphaa = 0.1;
    hold on
    d2 = [all_d, flip(all_d)];
    fill(d2,[prctile(time_count(2,:,:),prc,3), flip(prctile(time_count(2,:,:),100-prc,3))], ...
         [0, 0.447, 0.741],'LineStyle','none');
    alpha(alphaa)
    xlim([10,100]);
end
remove_border()
cd figs/
saveas(gcf, 'time', 'svg')
cd ..

figure
font_size = 18;
h2 = plot(all_d, mean(iter_count(2,:,:),3)', 'linewidth', 4);
grid on
set(gca, 'FontSize', 12);
ylabel('$\#$ of iterations','FontSize', font_size, 'Interpreter', 'latex');
xlabel('Dimensions $(d)$','FontSize', font_size, 'Interpreter', 'latex');
if shaded
    prc = 0;
    alphaa = 0.1;
    hold on
    d2 = [all_d, flip(all_d)];
    fill(d2,[prctile(iter_count(2,:,:),prc,3), flip(prctile(iter_count(2,:,:),100-prc,3))], ...
         [0, 0.447, 0.741],'LineStyle','none');
    alpha(alphaa)
    xlim([10,100]);
end
remove_border()
cd figs/
saveas(gcf, 'iteration', 'svg')
cd ..
