clear all; close all; clc;
%%
rng(0,'twister');

% Number of periods in the model
T = 100;

rho_vect = [0, 5e-3:5e-3:0.06];
c_vect = [0, 5e-3:5e-3:0.06];
delta = 1e-4;

% generate the seed pool of 1000 numbers
seedvalue = randperm(1000);

% set the number of runs
J = 200;


% Initialize the state of the inverted pendulum
initial_state = [0; 0; 0.1; 0.5];
[input.system, input.s, input.belief] = initialize(initial_state);

Wass_true_state = {};
Wass_estimated_state = {};
parfor i = 1:length(rho_vect)
    rho = rho_vect(i);
    display(['Wass iteration: ', num2str(i)]);
    
    for j = 1:J
        [ts, es] = EpisodeKalman( input, T, rho, delta, seedvalue(j), 'w');
        Wass_true_state{i, j} = ts;
        Wass_estimated_state{i, j} = es;
    end
end
labBarrier

kl_true_state = {};
kl_estimated_state = {};
parfor i = 1:length(c_vect)
    rho = c_vect(i);
    display(['KL iteration: ', num2str(i)]);
    for j = 1:J        
        [ts, es] = EpisodeKalman( input, T, rho, delta, seedvalue(j), 'kl');
        kl_true_state{i, j} = ts;
        kl_estimated_state{i, j} = es;
    end
end
labBarrier

% v_true_state = {};
% v_estimated_state = {};
% disp('Vanilla run');
% parfor j = 1:J    
%     [ts, es] = EpisodeKalman(input, T, 0, delta, seedvalue(j), 'vanilla');
%     v_true_state{j} = ts;
%     v_estimated_state{j} = es;
% end
% labBarrier

disp('Simulation finished!');
%disp('Saving data!');
%filename = 'Result.mat'
%save(filename);


%% Analyze the results from Kalman filtering
% Plot average MSE at time 51 (equivalent to t=50 in the paper)
start_time = 51;
end_time = 51;
k = 3;          % index of the state variable to track: [s, sdot, theta, thetadot]
J = 200;
w_mse = zeros(length(rho_vect), end_time);
kl_mse = zeros(length(c_vect), end_time);
v_mse = zeros(end_time, 1);

for i = 1:length(rho_vect)
    for j = 1:J
        ts = Wass_true_state{i, j};
        es = Wass_estimated_state{i, j};
        for t = start_time:end_time
            w_mse(i, t) = w_mse(i, t) + (ts(k, t) - es(k,t))^2;
        end
    end
end

for i = 1:length(c_vect)
    for j = 1:J
        ts = kl_true_state{i, j};
        es = kl_estimated_state{i, j};
        for t = start_time:end_time
            kl_mse(i, t) = kl_mse(i, t) + (ts(k, t) - es(k,t))^2;
        end
    end
end

figure;
hold on;

semilogx(c_vect, 10*log10(nanmean(kl_mse(:,start_time:end_time),2)), '-*', 'LineWidth', 2)
semilogx(rho_vect, 10*log10(nanmean(w_mse(:,start_time:end_time),2)), '-d', 'LineWidth', 2)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
legend({'Kullback-Leibler KF', 'Wasserstein KF'}, 'Interpreter', 'latex', 'Location', 'NW', 'FontSize', 16);
xlabel('$\rho$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Empirical Error (dB)', 'Interpreter', 'latex', 'FontSize', 20);