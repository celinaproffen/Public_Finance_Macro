% define parameters
% ------------------------------------------------------------------------

global alpha beta tau lambda T N eta varrho zeta y
alpha  = 0.3;
beta   = 0.99^40;
tau    = 0;
lambda = 0.5;
T      = 100;
N      = 1000; % number of simulations 

% set seed
rng(2)

% define shocks (in function??)
% ------------------------------------------------------------------------

% define shocks (temp!!)
eta    = [0.26, 0.407, 0.554, 0.701, 0.848, 0.995, ...
            1.142, 1.289, 1.436, 1.583, 1.73];
varrho = [0.54, 1.46];
zeta   = [0.87, 1.13];
y      = [10, 20];  % income outcome:[y_low, y_high]

% simulation of capital accumulation
% ------------------------------------------------------------------------

% aggregate shocks
[shock, ssl] = shockaggr(eta, varrho, zeta);

% find optimal saving rate given shocks
[sav] = saving(shock, ssl);

% create the shock history for y and zeta
[shockhistz, shockhisty] = shockhistfunc(zeta, y);

% compute starting (steady state) value of k (k_0)
logk0 = (log((1 - alpha)/(1 + lambda)) + log(sav))/(1-alpha);

% compute capital accumulation history
khist = khistfunc(logk0, sav, shockhisty, shockhistz);

% plot
plot(khist)



% ---------------------------------------------------------------------- %
%                             Functions                                  %
% ---------------------------------------------------------------------- %

% shock grid construction function
% ------------------------------------------------------------------------
function [shock, ssl] = shockaggr(eta, varrho, zeta)
    shock = combvec(eta, varrho, zeta)';
    ssl = length(shock); % state space length
end


% Saving function
% -------------------------------------------------------------------------
function sav = saving(shock, ssl)
    global beta alpha lambda tau

    phi = zeros(ssl, 1);
    for n = 1:ssl
        phi(n, 1) = 1 / (1 + (1 - alpha) / (alpha * (1 + lambda) * ... 
            shock(n, 2)) * (lambda * shock(n, 1)) + ...
            tau * (1 + lambda * (1 - shock(n, 2))));
    end
    
    Phi = sum(phi)/ssl;
    sav = (beta * Phi) / (1 + beta * Phi);
end

% create a random sequence of shocks
% ------------------------------------------------------------------------
function [shockhistz, shockhisty] = shockhistfunc(zeta, y)
    global T
    
    shockhistz = zeros(T, 1);
    shockhisty = zeros(T, 1);

    for i = 1:T
        shockhistz(i) = zeta(randperm(length(zeta), 1));
        shockhisty(i) = y(randperm(length(y), 1));
    end
end

% compute capital accumulation history given starting values and shock hist. 
% ------------------------------------------------------------------------
function khist = khistfunc(logk0, sav, shockhisty, shockhistz)

    global T

    khist = zeros(T, 1);
    
    % define starting value
    khist(1) = kapaccfunc(logk0, sav, shockhisty(1), shockhistz(1));

    % computing k_t t = 1:50000
    for i = 2:T
        khist(i) = kapaccfunc(khist(i), sav, shockhisty(i), shockhistz(i));
    end
    
    % logcapital law of motion 
    % --------------------------------------------------------------------
    function logk_1 = kapaccfunc(logk_0, sav_rate, y_1, zeta_1)
        global alpha lambda tau

        logk_1 = log((1 - alpha)/(1 + lambda)) + log(sav_rate) + ...
            log(y_1 - tau) + log(zeta_1) + alpha * logk_0;
    end
end
