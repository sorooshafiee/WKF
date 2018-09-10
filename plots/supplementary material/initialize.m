function [ system, s, belief ] = initialize(initial_state)
    system = initializeSystem;
    s = initializeState(initial_state);
    belief = initializeBelief(s);
    
    cost = initializeCost;
    [system.K, S, E] = dlqr(system.A, system.B, cost.Q, cost.R);
end

  
function [cost] = initializeCost()
    cost.Q = eye(4);
    cost.Q(1, 1) = 1e-2; % cost for x
    cost.Q(2, 2) = 1e-2; % cost for xdot
    cost.Q(3, 3) = 1e1; % cost for theta
    cost.Q(4, 4) = 1e-1; % cost for theta dot
    
    cost.R = 1e-5*eye(1);
    
    cost.Qf = eye(4);
    cost.Qf(1, 1) = 1; % cost for x
    cost.Qf(2, 2) = 1; % cost for xdot
    cost.Qf(3, 3) = 1e1; % cost for theta
    cost.Qf(4, 4) = 1e0; % cost for theta dot
end
    
  
function [ system ] = initializeSystem()
    % Generate the system for the cart pole
    system.g               = 9.8;     % Gravitational constant
    system.Mass_Cart       = 1;      % Mass of the cart is assumed to be 1Kg
    system.Mass_Pole       = 0.1;      % Mass of the pole is assumed to be 0.1Kg
    system.Total_Mass      = system.Mass_Cart + system.Mass_Pole;
    system.Length          = 0.3;      % Distance from pivot to center of mass of the rod
    system.PoleMass_Length = system.Mass_Pole * system.Length;
    system.Tau             = 0.1;     % Time interval for updating the values
    
    system.b               = 0.01;      % Friction coefficient
    system.I               = system.Mass_Pole*system.Length*system.Length/3;
    
    system.n               = 4;         % number of states
    
    denom = ((system.Mass_Pole + system.Mass_Cart)*(system.Mass_Pole*system.Length^2 + system.I) - system.Mass_Pole^2*system.Length^2);
    temp1 = - system.b*(system.Mass_Pole*system.Length^2 + system.I)/denom;
    temp2 = - system.Mass_Pole^2*system.Length^2*system.g/denom;
    temp3 = - system.b*system.Mass_Pole*system.Length/(-denom);
    temp4 = - (system.Mass_Pole + system.Mass_Cart)*system.Mass_Pole*system.g*system.Length/(-denom);
    cont_A = [0, 1, 0, 0; 0, temp1, temp2, 0; 0, 0, 0, 1; 0, temp3, temp4, 0]; 
    
    temp5 = (system.Mass_Pole*system.Length^2 + system.I)/denom;
    temp6 = (system.Mass_Pole*system.Length)/(-denom);
    cont_B = [0; temp5; 0; temp6];
    
    % set the C and D matrix for the observation
    %cont_C = eye(4); system.m = 4; cont_D = zeros(4, 1); % observe all 4 states
    cont_C = zeros(2, 4); cont_C(1, 1) = 1; cont_C(2, 3) = 1; system.m = 2; cont_D = zeros(2, 1); % observe only x and theta
    
    
    %states = {'x' 'x_dot' 'theta' 'theta_dot'};
    %inputs = {'u'};
    %outputs = {'x'; 'theta'};
    %sys_ss = ss(cont_A, cont_B, cont_C, cont_D,'statename',states,'inputname',inputs,'outputname',outputs);
    
    % Create the continuous time system
    % s_dot = cont_A*s + cont_B*u;
    % y = cont_C*s + cont_D*u;
    sys_ss = ss(cont_A, cont_B, cont_C, cont_D);
    
    % Create the discrete time system with A, B, C, D
    discrete = c2d(sys_ss,system.Tau,'zoh');
    
    system.A = discrete.A;
    system.B = discrete.B;
    system.C = discrete.C;
    
    % Initialize matrix W and V as in paper
    system.W = [0.1*eye(system.n), zeros(system.n, system.m)];         
    system.V = [zeros(system.m, system.n), 0.1*eye(system.m)];        
    
    
    system.true_W = system.W + 0.02*randn(size(system.W));
    system.true_W(:, end-1:end) = 0;
    system.true_V = system.V + 0.02*randn(size(system.V));
    system.true_V(1:system.m, 1:system.n) = 0;
    
    system.d = system.m + system.n;
end


function [s] = initializeState(state)
% x,x_dot,theta,theta_dot
    s = state;
end

function [belief] = initializeBelief(s)
    % mean
    belief.m = zeros(size(s));
    
    % covariance
    belief.cov = 0.01*eye(4);
end