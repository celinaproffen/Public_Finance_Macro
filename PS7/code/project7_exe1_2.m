%% Proposition 3 from Harenberg and Ludwig (2015)
%% Equilibrium dynamics of econ. 
% Model with:
% - aggr. productivity shocks
% - shocks to user costs of capital
% - 100% depreciation

%% Inititiate the code
% ------------------------------------------------------------------------
% define parameters
global allpha betta ttau llambda sav_ss log_k_ss k_ss zet zetta varrho etta 
allpha  = 0.3;
betta   = 0.99^40;
ttau    = 0;
llambda = 0.5;

% define length of simulations
T      = 100; % Alternatively change to 50,000 - should work fine! :)

% set seed for iid shocks
rng(2)

% Aggregate state: zet can be low or high
% shocks: zetta varrho etta 
% zetta (wages; low or high), varrho (rates of return; low or high) 
% and etta (Gaussian quadrature methos, 11 nodes); 
% all normally distributed with mean 1 and SD as in PS

zet = [0, 1];
zetta   = [0.87, 1.13];
varrho = [0.54, 1.46];
etta = [-4.0000, -3.2192, -2.6497, -1.8838, -0.9782, ...
  -0.4500, 0.5006, 1.3807, 2.1251, 2.6785, 3.0000];

%% Initiate Simulation

% Equilibrium Savings rate (equation 10 and 11 from Hatenberg & Ludwig)
Phhi = (1+ ((1-allpha)*(llambda + ttau*(1)))/(allpha*(1+llambda)))^(-1);
betta_Phhi = betta * Phhi;
sav_ss = betta_Phhi/(1+betta_Phhi);
log_k_ss = (1/(1-allpha))*(log((1-allpha)/(1+llambda))+log(sav_ss)+log(1-ttau));
k_ss = exp(log_k_ss);
%close all

[k_history, sav_history] = stoch_simulation(T)

%% Simulation of first order difference equation

% Create a vector of shocks of length T; whenever the economy is in a boom 
% or a recession (i.e. zet_low =0, zet_high = 1), zetta and varrho also
% need to be in that sstate

function [k_hist, sav_hist] = stoch_simulation(TT)
    %close all
    global allpha betta ttau llambda sav_ss log_k_ss k_ss zet zetta varrho etta 
    
    % initiate the whole thing
    log_k_hist = zeros(TT, 1);
    %k_hist = zeros(TT, 1);
    sav_hist = zeros(TT, 1);
    % Create shock series: column 1: zetta, column 2: varrho, column 3: etta
    shock_series = zeros(TT, 3); 
    
    % Draw random shock realizations from a combination
    shock_help = combvec(zet, etta)'; % find all combinations of aggregate and etta shocks 
    shock_combinations = length(shock_help); % state space length
    % Example, randi(10,[3,4]) returns a 3-by-4 array of pseudorandom integers between 1 and 10.
    shock_index_help = randi(shock_combinations,[TT,1]);
    
    for i = 1:TT
        a = shock_index_help(i,1)
        
        % Assign etta its value easily
        shock_series(i,3) = shock_help(a,2);
        
        if shock_help(a,1) == 0;
            shock_series(i,1) = zetta(1);
            shock_series(i,2) = varrho(1);
        elseif shock_help(a,1) == 1;
            shock_series(i,1) = zetta(2);
            shock_series(i,2) = varrho(2)
        else;
            disp('we seem to have a problem with the random number generator');
        end;
    end
    
    log_k_hist(1) = log_k_ss
    for j = 2:TT
        log_k_hist(j) = (log((1-allpha)/(1+llambda))+log(sav_ss)+ ...
            log(1-ttau)) + log(shock_series(j-1,1)) + allpha*log_k_hist(j-1);
    end
    
    k_hist = exp(log_k_hist)
    % plot:
    plot(k_hist)
end

%% 
