function [xhat, V, G, S] = WKF(sys, rho, Y_T, x_0, V_0, opts)
%   Wasserstein Kalman filter - MATLAB implementation of Algorithm 3
%
%   Syntax: [xhat, V, G, S] = WKF(sys, rho, Y_T, x0, V0, opts)
%   WKF() computes the WKF state estimate and optimal Nash equilibrium
%
%   Consider the linear dynamic system
%           x_t = A_t x_{t-1} + B_t v_t
%           y_t = C_t x_t + D_t v_t
%   where x_t is the unobserved state, y_t is the observed output, and
%   v_t is the stochastic driving noise of the system.
%
%   WKF has the following input arguments:
%   sys:  A structure contains the linear system model
%           sys.A:  Matrix A_t in state-space representation
%                   sys.A has to be a matrix with size (n * n) or
%                   a vector cell with length T with elments whose size is (n * n)
%           sys.B:  Matrix B_t in state-space representation
%                   sys.B has to be a matrix with size (n * d) or
%                   a vector cell with length T with elments whose size is (n * d)
%           sys.C:  Matrix C_t in state-space representation
%                   sys.C has to be a matrix with size (m * n) or
%                   a vector cell with length T with elments whose size is (m * n)
%           sys.D:  Matrix D_t in state-space representation
%                   sys.D has to be a matrix with size (m * d) or
%                   a vector cell with length T with elments whose size is (m * d)
%   rho:  A scaler or a vector with size T whose t'th element is the Wasserstein ambiguity size at time t
%   Y_T:  A matrix with size (m * T) whose t'th column is the observed output at time t
%   x0:   Initial point estimate
%   V0:   Initial covariance estimate
%   opts: Structure contains the parameters for Frank-Wolfe algorithm
%
%   WKF returns the following outputs:
%   xhat: A matrix with size (n * T) whose t'th column is the estimated state at time t
%   V:    A matrix with size (n * n * T) whose t'th block (V(:,:,t)) is the estimated
%         conditional covariance at time t
%   G:    A matrix with size (n * m * T) whose t'th block (G(:,:,t)) is the WKF gain at time t 
%   S:    A matrix with size (d * d * T) whose t'th block (S(:,:,t)) is the least favorable
%         covariance matrix at time t

    if nargin < 6
        opts = [];
    end
    are_cells = false;
    if iscell(sys.A) && iscell(sys.B) && iscell(sys.C) && iscell(sys.D) 
        T = min(length(sys.A), size(Y_T, 2));
        n = size(sys.A{1}, 1);
        m = size(sys.C{1}, 1);
        d = n + m;
        are_cells = true;
    else
        if iscell(sys.A) || iscell(sys.B) || iscell(sys.C) || iscell(sys.D)
            error('The system parameters should be either cell or matrix.')
        end
        T = size(Y_T, 2);
        n = size(sys.A, 1);
        m = size(sys.C, 1);
        d = n + m;
    end

    x_prev = x_0;
    V_prev = V_0;
    xhat = zeros(n, T);
    V = zeros(n, n, T);
    G = zeros(n, m, T);
    S = zeros(d, d, T);
    for t = 1 : T

        if are_cells
            [mu_t, Sigma_t] = predict(x_prev, V_prev, sys.A{t}, sys.B{t}, sys.C{t}, sys.D{t});
        else
            [mu_t, Sigma_t] = predict(x_prev, V_prev, sys.A, sys.B, sys.C, sys.D);
        end

        try
            rho_t = rho(t);
        catch
            rho_t = rho;
        end

        [xhat(:,t), V(:,:,t), G(:,:,t), S(:,:,t)] = update(mu_t, Sigma_t, rho_t, Y_T(:,t), length(x_0), opts);

        x_prev = xhat(:,t);
        V_prev = V(:,:,t);

    end
end

function [mu_t, Sigma_t] = predict(x_prev, V_prev, A_t, B_t, C_t, D_t)
    A_aug = [A_t; C_t * A_t];
    B_aug = [B_t; C_t * B_t + D_t];
    mu_t = A_aug * x_prev;
    Sigma_t = A_aug * V_prev * A_aug' + B_aug * B_aug';
end

function [xhat_t, V_t, G_t, S_t] = update(mu_t, Sigma_t, rho_t, y_t, x_dim, opts)
    [phi_star, Q_star] = Frank_Wolfe(mu_t, Sigma_t, rho_t, x_dim, opts);
    G_t = phi_star.G;
    S_t = Q_star.Sigma;
    V_t = S_t(1:x_dim, 1:x_dim) - G_t * S_t(x_dim+1:end, 1:x_dim);
    xhat_t = G_t*(y_t - mu_t(x_dim+1:end,1)) + mu_t(1:x_dim,1);
end