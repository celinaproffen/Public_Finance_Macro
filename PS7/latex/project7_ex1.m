global alpha beta tau lambda T N s phi Phi eta varrho zeta shock ssl
alpha = 0.3;
beta = 0.99^40;
tau = 0;
lambda = 0.5;
T = 50000;
N = 1000; % number of simulations
s = zeros(1, T);
phi = zeros(ssl, 1);
Phi = zeros(1, 1);

% solving 
%-------------------------------------------------------------------------

% find optimal saving given shocks
[Phi] = phi_func(beta, alpha, tau, lambda, shock, ssl);
[s] = saving(beta, Phi);
logs = log(s);

% computing starting value of k (k_0)
k0 = (log((1 - alpha)/(1 + lambda)) + logs)/(1-alpha);

% create the shock history for t = 1:50,000 and z and zeta

    % draft 50000 value of z
        
        % define
        zl = 10;
        zh = 20;
        
    % draft 50000 value of zeta

% computing k_t t = 1:50000

% qes

% function


% set up grid
%-------------------------------------------------------------------------

% set up function

eta = zeros(11, 1);
varrho = zeros(2, 1);
zeta = zeros(2, 1);
shock = zeros(44, 3);  

%...
% output
eta = [0.26, 0.407, 0.554, 0.701, 0.848, 0.995, 1.142, 1.289, 1.436, 1.583, 1.73];
varrho = [0.54, 1.46];
zeta = [0.87, 1.13];

% create shocks grid
shock = combvec(eta, varrho, zeta)';
ssl = 44; % state space length

% Saving
%-------------------------------------------------------------------------
function s = saving(beta,Phi)
    s=(beta*Phi)/(1+beta*Phi)
end

 % function 
function Phi = phi_func(beta, alpha, tau, lambda, shock, ssl)
        for n = 1:ssl
            phi(n,1) = 1 / (1 + (1 - alpha)/(alpha * (1 + lambda) * shock(n, 2)) * (lambda * shock(n, 1)))
        end
        Phi = sum(phi)/ssl;
end
