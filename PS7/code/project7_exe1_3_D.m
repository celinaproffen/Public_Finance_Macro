%% Proposition 3 from Harenberg and Ludwig (2015)
%% Equilibrium dynamics of econ. 
% Model with:
% - aggr. productivity shocks
% - shocks to user costs of capital
% - 100% depreciation

% We want to evaluate welfare over different values of TAU
% We want to compute the consumption equivalent variation;
% For this We need 


% % Assignment Repeat the exercise for  = 0:1. Compute welfare in the two
% economies with tau = 0.0 and tau = 0.1 by taking an ex-ante perspective.
% To this purpose take the average (again discarding the first 500 simulations) of realized cohort utilities 
% Compute the consumption equivalent variation from the ex-ante perspective and interpret your findings.

%% Inititiate the code
% ------------------------------------------------------------------------
% define parameters
global allpha betta llambda sav_ss log_k_ss k_ss zet zetta varrho etta tech
allpha  = 0.3;
betta   = 0.99^40;
ttau_orig = 0;
ttau_altern = 0.1;
llambda = 0.5;
tech = 1; % technology level (not needed here yet, as it doesn't grow)

% define length of simulations
T      = 50000; % later 50,000
M      = 12; % number of iterations 
N      = 1; % number of agents

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
  -0.4500, 0.9506, 1.8307, 2.5751, 3.1285, 3.0000];

%% Initiate Simulation


%close all

% We want these functions to spit out the final Values under these regimes
% Hence I added a simulation which happens at the end or upon convergence
% and calculates the expectd utilities; 
% We should include the lifetime values of the households that are
% infinitely lived
V_orig = KrusselSmith(M,T,ttau_orig)
V_altern = KrusselSmith(M,T,ttau_altern)
eqiv = EQIV(V_altern,V_orig)

%% Simulation of first order difference equation (Exercise 1.2)

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
            shock_series(i,2) = varrho(2);
        else;
            disp('we seem to have a problem with the random number generator');
        end;
    end
    
    log_k_hist(1) = log_k_ss
    for j = 2:TT
        log_k_hist(j) = (log((1-allpha)/(1+llambda))+log(sav_ss)+ ...
            log(1-ttau)) + log(shock_series(j-1,1)) + allpha*log_k_hist(j-1);
    end
    
    k_hist = exp(log_k_hist);
    % plot:
    plot(k_hist)
end

%% Simple Krusell-Smith Algorithm (Exercise 1.3)
% --> now expanded to account for different tax policies

function expec_utility = KrusselSmith(M,T,tauu)
global allpha betta llambda sav_ss log_k_ss k_ss zet zetta varrho etta tech

ttau = tauu;
M = M;
T = T;

% Equilibrium Savings rate (equation 10 and 11 from Hatenberg & Ludwig)
Phhi = (1+ ((1-allpha)*(llambda + ttau*(1)))/(allpha*(1+llambda)))^(-1);
betta_Phhi = betta * Phhi;
sav_ss = betta_Phhi/(1+betta_Phhi);
log_k_ss = (1/(1-allpha))*(log((1-allpha)/(1+llambda))+log(sav_ss)+log(1-ttau));
k_ss = exp(log_k_ss);

% Compute the different values for psi_zero and psi_one given each state
% psi init has two rows: pssi_zero_init and pssi_one_init
pssi_zero_init = zeros(1,2);
pssi_one_init = zeros(1,2);

% Assuming that the expectations of future shocks equals zero, these are
% the theoretical moments:
    for i = 1:2
        pssi_zero_init(1,i) = (log((1-allpha)/(1+llambda))+log(sav_ss)+log(1-ttau)) + log(zetta(i));
        pssi_one_init(1,i)= allpha;
    end

psi_init = zeros(2,2);
psi_init(1,:) = pssi_zero_init;
psi_init(2,:) = pssi_one_init;

%% Iterate over the values of psi 
% (first column low state; second column high state)

% Adapt if needed:
gridsize_shock_help = size(zet,[2]); % number of aggregate shocks
gridsize_kapital_help = 5; % number of gridpoints for capital
ommega = 0.95;
tol = eps;

% Create grid for possible capital stock values
minn = max(0.5*k_ss,sqrt(eps));
maxx = 1.5*k_ss;
grid_kapital = linspace(minn,maxx,gridsize_kapital_help);

% store saving rates on a grid with two dimensions (zet and k_grid)
grid_sav = zeros(gridsize_shock_help,gridsize_kapital_help);
    for i = 1:gridsize_kapital_help
        grid_sav(:,i)= zet';
    end

% grid of capital over iterations; iniatialize it with capital grid..
grid_kapital_dynamics = zeros(T,gridsize_shock_help,gridsize_kapital_help);
    for j = 1:T
        for i = 1:gridsize_kapital_help
            grid_kapital_dynamics(j,:,i)= zet';
        end
    end
    for i = 1: gridsize_shock_help
        grid_kapital_dynamics(1,i,:) = log(grid_kapital);
    end

% Draw random series of shocks: columns represent zetta, varrho, etta
[aggregate_shocks, zettavarrhoetta_shocks] = create_shocks(T);

% Create a grid to store the utilities of each cohort given shocks, initial capital etc  
U = zeros(T,gridsize_shock_help,gridsize_kapital_help);

% -------------------------------------------------------------------------
% starting value psi_init and capital_0
psi_new = psi_init;

    for m = 2:M
        
        psi_old = psi_new;
                     
        % Household Model from Harenberg and Ludwig (2015), equations 1-3
        % Assumption: u(c) = log(c)
        
        % find a saving rate associated with this
        % store the solution for the savings rate on a grid with two
        % dimensions: firstly, aggregate state (boom, recession) and
        % secondly, initial level of capital k_t
        
        % I need to define a function for the Euler equation. The solution
        % when this Euler equation is zero, is how the savings rate is defined
        % --> See Project description submitted as Latex file
        
        for i = 1:gridsize_shock_help
            
            for k = 1:gridsize_kapital_help
            
                % Define shock and initial capital you are looping over
                zetta_0 = zetta(i);
                kap_0 = grid_kapital(k);
                
                % Feed fzero the elements it needs to solve for the savings
                % rate
                myfun = @(sav,zetta_0,kap_0,ttau) log(betta) + log(sav^(allpha-1)*(1-ttau)^(allpha-1)*allpha*(1-allpha)^(allpha-1)*(1+llambda)^(1-allpha)*kap_0^(allpha*(allpha-1))*zetta_0^(allpha-1)) + ...
                        log((1-sav)*(1-ttau)*(1-allpha)*kap_0^(allpha)*zetta_0) - log(sav^(allpha)) - ...
                        log((1-ttau)^(allpha)*allpha*(1-allpha)^(allpha-1)*(1+llambda)^(1-allpha)*kap_0^(allpha*(1-allpha))*1^(allpha-1) + ...
                        (1-ttau)^(1+allpha)*llambda*(1+llambda)^(-allpha)*(1-allpha)^(1+allpha)*kap_0^(allpha*allpha)*1^(allpha) + ...
                        ttau*(1-ttau)^(allpha)*(1-llambda)*(1+llambda)^(-allpha)*(1-allpha)^(1+allpha)*kap_0^(allpha*(allpha-1))*1^(allpha))   ;  % parameterized function
                fun = @(sav) myfun(sav,zetta_0,kap_0,ttau);    % function of savings rate alone
                x0 = [0.01 0.99]; % initial interval; the savings rate should by definition be between zero and one
                %options = optimset('Display','iter'); % show iterations
                               
                % solve for the savings rate
                %  !!! Array indices must be positive integers or logical values.
                [sav_help] = fzero(fun,0.5); %fzero(fun,0.5)
                
                % Save the value of savings rate for aggregate state and initial capital             
                grid_sav(i,k) = sav_help;
                
                % Current period 1 utility
                U(1,i,k) = log((1- grid_sav(i,k))*(1-allpha)*kap_0^(allpha)*zetta_0);
                
                
                % Below:
                % A) Iterate forward for capital given the savings rate and the realization of shocks
                % B) Calculate future realized contemporaneous utilities  
                
                for t = 2:T
                    
                    grid_kapital_dynamics(t,i,k) = (log((1-allpha)/(1+llambda)) + log(grid_sav(i,k))+ ...
                               log(1-ttau)) + log(zettavarrhoetta_shocks(t-1,1)) + allpha*grid_kapital_dynamics(t-1,i,k) ;
                    
                    U(t,i,k) = log((1- grid_sav(i,k))*(1-allpha)*exp(grid_kapital_dynamics(t,i,k))^(allpha)*zettavarrhoetta_shocks(t,1));
                        
                end % end T loop
                
            end % end K loop
            
            % Delete first 500 observations for each shock state and initial capital
            dellete = 500;
            grid_kapital_dynamicshelp = grid_kapital_dynamics(dellete+1:end,:,:); % Extracts rows (501) until the bottom.
            
            % OBJECTIVE 1: FIND PSI
            X = squeeze(grid_kapital_dynamicshelp(1,i,:)); 
            y = squeeze(grid_kapital_dynamicshelp(2,i,:)); 

            % For each shock state and level of starting capital define regressor and regressand  
            for tt = 3:(T-dellete)
                X = [X;squeeze(grid_kapital_dynamicshelp(tt-1,i,:))];
                y = [y;squeeze(grid_kapital_dynamicshelp(tt,i,:))];
            end % end of building regressor and regressand
            
            % obtain new guess of psi's (psi_one and psi_zero for each state of the economy)
            % i.e. obtain psi_new (note: fitlm includes a constant automatically ;)
            lm = fitlm(X,y);
            psi_new(1,i) = lm.Coefficients.Estimate(1); % psi_zero coefficient
            psi_new(2,i) = lm.Coefficients.Estimate(2); % psi_one coefficient
            
            clear X y
                        
        end % end i loop
        
        % OBJECTIVE 2: CALCULATE EXPECTED UTILITY
        % Average over time
        % For each point in time: 
        % Average over initial level of capital
        % Average over shocks
        expec_utility_help = U(dellete+1:end,:,:); % Extracts rows (501) until the bottom.
        expec_utility_help2 = mean(expec_utility_help,[2 3]);
        expec_utility_help3 = squeeze(expec_utility_help2(:,1,1));
        % sum and discount all this
        betta_help = ones(T-dellete,1);
            for ttt = 1:T-dellete
                betta_help(ttt) = betta^(ttt);
            end
        expec_utility = (expec_utility_help3')*betta_help;
        
        % compare new guess of psi to old one
        diff_psi = psi_new - psi_old;
        
            if (abs(diff_psi)<tol),
                disp('convergence');
                disp(['final PSI:', num2str(psi_new)]);
                break;
             else,
                 % update guesses of psi and capital
                 psi_new = ommega*psi_old +(1-ommega)*psi_new;
                 disp(['iteration #', num2str(m)]);
            end % end updating psi

    end % end M loop
    
    % Let the function display the initial and final version of PSI
    disp(['Initial theoretical value of PSI:']);
    disp(psi_init)
    disp(['Final regressed value of PSI:']);
    disp(psi_new)
   
end  % end function


%% Function to create random draws for the shock series of zetta, varrho and etta
function[zet_series,shock_series] = create_shocks(iterations)

    global allpha betta ttau llambda sav_ss log_k_ss k_ss zet zetta varrho etta 

    shock_series = zeros(iterations, 3); 
    
% Draw random shock realizations from a combination
    shock_help = combvec(zet, etta)'; % find all combinations of aggregate and etta shocks 
    shock_combinations = length(shock_help); % state space length
    % Example, randi(10,[3,4]) returns a 3-by-4 array of pseudorandom integers between 1 and 10.
    shock_index_help = randi(shock_combinations,[iterations,1]);
    zet_series = zeros(iterations,1);
    
    for i = 1:iterations
        a = shock_index_help(i,1);
        
        % Assign etta its value easily
        shock_series(i,3) = shock_help(a,2);
        
        if shock_help(a,1) == 0;
            zet_series(i) = 0;
            shock_series(i,1) = zetta(1);
            shock_series(i,2) = varrho(1);
        elseif shock_help(a,1) == 1;
            zet_series(i) = 1;
            shock_series(i,1) = zetta(2);
            shock_series(i,2) = varrho(2);
        else;
            disp('we seem to have a problem with the random number generator');
        end
    end
end

%% Equivalent variation
function eqiv=EQIV(V_altern,V_orig)
    %global tetta
    %eqiv = (V/V_bar)^(1/(1-tetta))-1;
    
end

