function [ true_state, estimated_state] = EpisodeKalman(input, maxsteps, rho, delta, seedvalue, type)

    rng(seedvalue, 'twister');

    system = input.system;
    s = input.s;
    belief = input.belief;

    % a matrix to keep track of all the states over time
    true_state = nan(system.n, maxsteps+1);
    estimated_state = nan(system.n, maxsteps+1);

    for i=1:maxsteps  
        %display(['Step:', num2str(i)]);
        true_state(:, i) = s;
        estimated_state(:, i) = belief.m;

        % Find the optimal control as a linear feedback of the state estimate
        % notice that by the optimal action has a negative sign 
        % For a reason of the negative sign, type: help dlqr
        action = - system.K*belief.m;

        % do the selected action and get the true state    
        sp  = TrueStateEvolution( system, action , s );    

        % add some random noise to the true state
        sp = sp + 0.5*system.true_W*randn(system.d, 1);

        
        % Do Kalman filtering
        belief = KalmanPredict(system, belief, action);

        % Find the worst-case covariance matrix
        switch type
            case 'vanilla'
                % Vanilla Kalman filter, or classical Kalman filter
                S = belief.predict_cov;
                G = S(1:system.n, system.n+1:end)/(S(system.n+1:end, system.n+1:end));
            case 'kl'
                [G, S] = tau_update(belief.predict_cov, rho, 0, system.n);
            case 'w'
                [phi_star, Q_star] = Frank_Wolfe(zeros(system.d,1), belief.predict_cov, rho, system.n);
                G = phi_star.G;
                S = Q_star.Sigma;
            otherwise
                display('Type of estimator unknown!');
        end
        % Generate an observation of the new true state
        noise = randn(system.d, 1);
        y = GetObservation(system, sp, noise);
        belief = KalmanCorrect(belief, G, S, system.n, y);

        %update the true state
        s = sp;
    end

    true_state(:, i+1) = s;
    estimated_state(:, i+1) = belief.m;
end

function [s] = TrueStateEvolution(system, force, s)
% True state evolution of the system

        % Parameters for simulation
        x          = s(1);
        x_dot      = s(2);
        theta      = s(3);
        theta_dot  = s(4);

        denom = ((system.Mass_Pole + system.Mass_Cart)*(system.Mass_Pole*system.Length^2 + system.I) - system.Mass_Pole^2*system.Length^2*cos(theta)^2);
        xacc = ((force - system.b*x_dot + system.Mass_Pole*system.Length*theta_dot*theta_dot*sin(theta))*(system.Mass_Pole*system.Length*system.Length + system.I) - system.Mass_Pole^2*system.Length^2*system.g*sin(theta))/denom;
        thetaacc = ((force - system.b*x_dot + system.Mass_Pole*system.Length*theta_dot*theta_dot*sin(theta))*(system.Mass_Pole*system.Length*cos(theta)) - (system.Mass_Pole + system.Mass_Cart)*system.Mass_Pole*system.g*system.Length*sin(theta))/(-denom);

        % Update the four state variables, using Euler's method.
        x         = x + system.Tau * x_dot;
        x_dot     = x_dot + system.Tau * xacc;
        theta     = theta + system.Tau * theta_dot;
        theta_dot = theta_dot + system.Tau*thetaacc;

        % return true future state
        s = [x; x_dot; theta; theta_dot];
end

function [ belief ] = KalmanPredict( system, belief, action )
    % Do the one-step ahead prediction step
    
    % In this version of the Kalman filter, the prediction is done over the
    % (x, y) space
    belief.predict_m    = [system.A*belief.m + system.B*action; system.C*(system.A*belief.m + system.B*action)];
    
    mat1 = [system.A; system.C*system.A];
    mat2 = [system.W; system.C*system.W + system.V];
    
    belief.predict_cov  = mat1*belief.cov*mat1' + mat2*mat2';

end


function [ belief ] = KalmanCorrect( belief, G, S, n, y )
    mean_x = belief.predict_m(1:n);
    mean_y = belief.predict_m(n+1:end);
    
    Sxx = S(1:n, 1:n);
    Syx = S(n+1:end, 1:n);
    
    
    % Update the belief about the current state
    belief.m = mean_x + G*(y - mean_y);
    belief.cov = Sxx - G * Syx;
end


function [ y ] = GetObservation( system, s, noise )
    y = system.C*s;
    y = y + system.true_V*noise;
end
