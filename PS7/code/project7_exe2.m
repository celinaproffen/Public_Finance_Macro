%% Proposition 3 from Harenberg and Ludwig (2015)
%% Equilibrium dynamics of econ. 
% Model with:
% - aggr. productivity shocks
% - shocks to user costs of capital
% - 100% depreciation

%% Inititiate the code
% ------------------------------------------------------------------------
% define parameters
global allpha betta ttau llambda sav_ss log_k_ss k_ss zet zetta varrho etta tech
allpha  = 0.3;
betta   = 0.99^40;
ttau    = 0;
llambda = 0.5;
tech = 1; % technology level (not needed here yet, as it doesn't grow)

% define length of simulations
T      = 50000; % later 50,000
M      = 100; % number of iterations 
N      = 1; % number of agents

% set seed for iid shocks
rng(2)

% Aggregate state: zet can be low or high
% shocks: zetta varrho etta 
% zetta (wages; low or high), varrho (rates of return; low or high) 
% and etta (Gaussian quadrature methos, 11 nodes); 
% all normally distributed with mean 1 and SD as in PS

zet = [0, 1];
zetta   = [0.97, 1.03];
varrho = [0.54, 1.46];
etta = [-4.0000, -3.2192, -2.6497, -1.8838, -0.9782, ...
  -0.4500, 0.9506, 1.8307, 2.5751, 3.1285, 3.0000];

%% Initiate Simulation

% Equilibrium Savings rate (equation 10 and 11 from Hatenberg & Ludwig)
Phhi = (1+ ((1-allpha)*(llambda + ttau*(1)))/(allpha*(1+llambda)))^(-1);
betta_Phhi = betta * Phhi;
sav_ss = betta_Phhi/(1+betta_Phhi);
log_k_ss = (1/(1-allpha))*(log((1-allpha)/(1+llambda))+log(sav_ss)+log(1-ttau));
k_ss = exp(log_k_ss);
%close all

%[k_history, sav_history] = stoch_simulation(N, T)
ComplexKrusselSmith(M,T)


%% Complex Krusell-Smith Algorithm (Exercise 1.4)

function ComplexKrusselSmith(M,T)
global allpha betta ttau llambda sav_ss log_k_ss k_ss zet zetta varrho etta tech

M = M;
T = T;

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

% Draw a random realization of shocks (necessary?)

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

% NEW: Draw random series of shocks DEPENDENT on initia state of the world
% columns of shock_series represent zetta, varrho, etta

% If we start from recession:
[low_aggregate_shocks, low_aggregate_shocks_probability, low_zettavarrhoetta_shocks] = create_Markov_shocks(T,0);
% If we start from boom:
[high_aggregate_shocks, high_aggregate_shocks_probability, high_zettavarrhoetta_shocks] = create_Markov_shocks(T,1);

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
        % when this Euler equation is zero, is how the savings rate is
        % define
        
        for i = 1:gridsize_shock_help
            
            for k = 1:gridsize_kapital_help
            
                % Define shock and initial capital you are looping over
                kap_0 = grid_kapital(k);
                zetta_0 = zetta(i);
                    if i ==1
                        zettavarrhoetta_shocks = low_zettavarrhoetta_shocks;
                    elseif i ==2
                        zettavarrhoetta_shocks = high_zettavarrhoetta_shocks;
                    else
                        disp('Again something got wrong with the shock indices..')
                    end
                                
                % Feed fzero the elements it needs to solve for the savings rate
                
                myfun = @(sav,zetta_0,kap_0) log(betta) + log(sav^(allpha-1)*(1-ttau)^(allpha-1)*allpha*(1-allpha)^(allpha-1)*(1+llambda)^(1-allpha)*kap_0^(allpha*(allpha-1))*zetta_0^(allpha-1)) + ...
                        log((1-sav)*(1-ttau)*(1-allpha)*kap_0^(allpha)*zetta_0) - log(sav^(allpha)) - ...
                        log((1-ttau)^(allpha)*allpha*(1-allpha)^(allpha-1)*(1+llambda)^(1-allpha)*kap_0^(allpha*(1-allpha))*zetta_0^(allpha-1) + ...
                        (1-ttau)^(1+allpha)*llambda*(1+llambda)^(-allpha)*(1-allpha)^(1+allpha)*kap_0^(allpha*allpha)*zetta_0^(allpha) + ...
                        ttau(1-ttau)^(allpha)*(1-llambda)*(1+llambda)^(-allpha)*(1-allpha)^(1+allpha)*kap_0^(allpha*(allpha-1))*zetta_0^(allpha-1))   ;  % parameterized function
                fun = @(sav) myfun(sav,zetta_0,kap_0);    % function of savings rate alone
                x0 = [0 1]; % initial interval; the savings rate should by definition be between zero and one
                %options = optimset('Display','iter'); % show iterations
                
                
                % solve for the savings rate
                [sav_help fval] = fzero(fun,0.5); %fzero(fun,x0)
                
                % Save the value of savings rate for aggregate state and initial capital             
                grid_sav(i,k) = sav_help;
                
                for t = 2:T
                    % Iterate forward for capital given the savings rate and the realization of shocks
                    
                    grid_kapital_dynamics(t,i,k) = (log((1-allpha)/(1+llambda)) + log(grid_sav(i,k))+ ...
                               log(1-ttau)) + log(zettavarrhoetta_shocks(t-1,1)) + allpha*grid_kapital_dynamics(t-1,i,k) ;
                        
                end % end T loop
                
            end % end K loop
            
            % Delete first 500 observations for each shock state and initial capital
            dellete = 500;
            grid_kapital_dynamicshelp = grid_kapital_dynamics(dellete+1:end,:,:); % Extracts rows (501) until the bottom.
            
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
    disp(['Final regressed value of PSI (left: recession in t= 0; right: boom in t = 0 ):']);
    disp(psi_new)

end  % end function


%% Function to create random draws for the shock series of zetta, varrho and etta
% This time it all depends on the initial state and the transition -->
% probabilities are not random; instead I added Markov-chain reliant
% changes using the following ideas from: 
% - https://de.mathworks.com/matlabcentral/answers/57961-simulating-a-markov-chain
% - https://de.mathworks.com/help/econ/dtmc.simulate.html

function[zet_series,zet_series_probability,shock_series] = create_Markov_shocks(iterations,start_state_agg)

    global allpha betta ttau llambda sav_ss log_k_ss k_ss zet zetta varrho etta 
    
    % Create Markov process
    P = [0.95 0.05; 0.05 0.95];
    mc = dtmc(P);
    rng(2); % For reproducibility
    numSteps = iterations;
    
    % initial state to create Markov distribution
    if start_state_agg == 0
        init_state = [1 0];
        starting_value = 1; 
    elseif start_state_agg == 1
        init_state = [0 1];
        starting_value = 2; 
    else 
        disp('Problem specifying initial state of Markov process');
    end
    
    shock_series = zeros(iterations, 3); 
    zet_series = zeros(iterations, 1);
    zet_series_probability = zeros(iterations,2);
    
    for tt = 1:iterations
        zet_series_probability(tt,:) = init_state*(P^(tt-1));
    end

    chain = zeros(iterations,1);
    chain(1)=starting_value;
    for tt=2:iterations
        this_step_distribution = P(chain(tt-1),:);
        cumulative_distribution = cumsum(this_step_distribution);
        r = rand();
        chain(tt) = find(cumulative_distribution>r,1);
    end
    
    zet_shock_series = chain;
    
    etta_length = length(etta'); % state space length
    shock_index_help = randi(etta_length,[iterations,1]);   
    for i = 1:iterations
        % Assign etta its value easily
        shock_series(i,3) = etta(shock_index_help(i,1));
    end
    
    
    for ttt = 1:iterations
        
        if zet_shock_series(tt) == 1;
            zet_series(ttt) = 0;
            shock_series(ttt,1) = zetta(1);
            shock_series(ttt,2) = varrho(1);
        elseif zet_shock_series(tt) == 2;
            zet_series(ttt) = 1;
            shock_series(ttt,1) = zetta(2);
            shock_series(ttt,2) = varrho(2);
        else;
            disp('we seem to have a problem with the random number generator');
        end
    end
end
