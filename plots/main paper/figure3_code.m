clc
clearvars

% Set seed number
rng(12345);

% Auxilary function
vec = @(x) x(:);
out_of_sample_loss = @(G, S) vec([eye(size(G,1)), -G; -G', G'*G])' * vec(S);

% Problem setting
n = 80;
m = 20;
d = n + m;
run_count = 10000;
opts.verbose = false;

% Initialization
Bayes = zeros(1, run_count);
DRO = zeros(1, run_count);
MMSE = zeros(1, run_count);

for r = 1 : run_count

    fprintf('Running Iteration %d\n', r);

    A = randn(d);
    [R, ~] = eig(A + A');
    lambda = 0.1 + 9.9*rand(d,1);
    Sigma = R * diag(lambda) * R';
    Sigma_half = R * diag(sqrt(lambda)) * R';

    A_star = randn(d);
    [R_star, ~] = eig(A_star + A_star');
    lambda_star = rand(d,1);
    Delta_star = R_star * diag(lambda_star) * R_star';
    Delta_star_half = R_star * diag(sqrt(lambda_star)) * R_star';

    Sigma_star = (Sigma_half + Delta_star_half)^2;

    G_star = Sigma_star(1:n, n+1:end)/(Sigma_star(n+1:end, n+1:end));
    MMSE(r) = out_of_sample_loss(G_star, Sigma_star);

    G_Bayes = Sigma(1:n, n+1:end)/(Sigma(n+1:end, n+1:end));
    Bayes(r) = out_of_sample_loss(G_Bayes, Sigma_star);

    phi = Frank_Wolfe(zeros(d,1), Sigma, sqrt(d), n, opts);
    DRO(r) = out_of_sample_loss(phi.G, Sigma_star);

end

%%
figure
font_size = 18;
h1 = histogram(Bayes-MMSE,'Normalization','pdf'); 
hold on; 
h2 = histogram(DRO-MMSE,'Normalization','pdf');
set(gca, 'FontSize', 12);
ylabel('Probability density','FontSize', font_size, 'Interpreter', 'latex');
xlabel('Relative mean square error','FontSize', font_size, 'Interpreter', 'latex');
leg1 = legend({'${\rm Bayes}$', '${\rm Wasserstein}$'}, 'Interpreter', 'latex');
set(leg1,'FontSize',16);
if d == 10
    xlim([-0.1,2.0])
elseif d == 50
    xlim([0.8,5])
elseif d == 100
    xlim([2.5,8])
end
remove_border()
cd figs/
saveas(gcf, 'd_' + string(d), 'svg')
cd ..
