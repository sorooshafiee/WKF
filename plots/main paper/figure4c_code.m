clc
clearvars

% Set seed number
rng(12345);

% Problem setting
opts.verbose = true;
opts.tol = -inf;
opts.iter_max = 1e4;
run_count = 100;
d = 100;

% Initialization
rel_gap = zeros(opts.iter_max, run_count);

for r = 1 : run_count

    fprintf('Running Iteration %d\n', r);

    rho = sqrt(d);
    n = fix(d/5);

    A = randn(d);
    [R, ~] = eig(A + A');
    lambda = 0.1 + 9.9*rand(d,1);
    Sigma = R * diag(lambda) * R';

    [~, ~, obj, res] = Frank_Wolfe(zeros(d,1), Sigma, rho, n, opts);
    rel_gap(:,r) = abs(res) ./ obj;

end

%%
figure
r = 1 : opts.iter_max;
shaded = true;
font_size = 18;
loglog(r, mean(rel_gap,2), 'linewidth', 4);
grid on
set(gca, 'FontSize', 12);
ylabel('Relative duality gap','FontSize', font_size, 'Interpreter', 'latex');
xlabel('Iteration','FontSize', font_size, 'Interpreter', 'latex');
if shaded
    prc = 0;
    alphaa = 0.1;
    hold on
    r2 = [r, flip(r)];
    fill(r2,[prctile(rel_gap, prc, 2)', flip(prctile(rel_gap, 100-prc, 2))'], ...
         [0, 0.447, 0.741], 'LineStyle', 'none');
    alpha(alphaa)
end
remove_border()
cd figs/
saveas(gcf, 'convergence', 'svg')
cd ..
