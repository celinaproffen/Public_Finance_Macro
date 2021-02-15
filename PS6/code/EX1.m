% exercize
function exercise

close all

global alpha betta tetta delta replrate tau lambda 
global retage maxage L R N_age sr 

rho = 0.05;  % 0.01
delta = 0.05;
alpha = 0.33; % 0.33; % 0.1;
ret = 0.07;
tetta = 1.0;
betta = 1/(1+rho);
lambda = 1.0;

maxage = 80;        % corresponds to real age 100 == nj
retage = 45;        % corresponds to actual age 65 == jr
replrate = 0.0;     % ?????

% iteration parameters
tol = 1e-4;          % convergence criteria
maxit = 100;         % max n. of iteration
df = 0.1;            % smoothing parameter

mr = readfile([],'MR.txt',3);
sr = 1.0-mr(21:end,1);

% compute cohort size
[N_age,L_age,R_age] = func_pop(sr,maxage,retage);

N = sum(N_age);      % entire pop 
L = sum(L_age);      % working age pop
R = sum(R_age);      % retired pop

% Parameters for stochastic component
opt_det = false;       
opt_ny = 2; 

% cash-on-hand  grid parameters
nx = 30;             % # of grid-points
curv = 3.0;          % curvature of grid
grdfac = 40;         % scaling factor of saving grid

% Plot supply and demand
%-------------------------------------------------------------------------

% backup starting interest rate
retold = ret;

% create interest rate grid
retvec = [-delta+0.01:0.001:0.1]';
n = length(retvec);

% create asset grid
assvec = zeros(n,1);

% for each interest rate point, compute assets (== capital) supplied
for i = 1:n,
    ret = retvec(i);
    [fval, assvec(i)] = func_olg(ret);
end;

% graph of K(r) = optimal demand (input choice):
kvec = (alpha./(retvec+delta)).^(1/(1-alpha)).*L;

% plot
figure; hold on;
plot(retvec,assvec,'b-','LineWidth',3);
plot(retvec,kvec,'r--','LineWidth',3);
legend('demand', 'supply');
title('supply and demand in OLG economy');
print -depsc equil.eps


% solution of OLG model:
%-------------------------------------------------------------------------

% (re)set starting interest rate
ret = retold;

% iterate until r* s.t. => fval < tol 
for it=1:maxit,
    
    % given ret, compute differene between demand and supply r
    [fval, ass, Y] = func_olg(ret);
    
    % if fval(ret) < tol, STOP 
    if (abs(fval)<tol),
        disp('convergence');
        break;
        
    % fval(ret) > tol
    else,
        
        % update ret using df as a smoothing parameter
        ret = ret - df * fval;
        disp(['iteration #', num2str(it)]);
        disp(['guess of rate of return: ', num2str(ret)]);
    end;
end;
if (it>=maxit),
    warning('no convergence');
end;
disp(['equilibrium interest rate: ', num2str(ret)]);
disp(['equilibrium contribution rate: ', num2str(tau)]);
disp(['equilibrium capital output ratio: ', num2str(ass/Y)]);

end


% compute the size of each cohort given min, max and survival rate sr
% ------------------------------------------------------------------------
function [N,L,R]=func_pop(sr,maxage,retage);

N = zeros(maxage,1);
L = zeros(maxage,1);
R = zeros(maxage,1);
N(1) = 100;
for age=2:maxage,
    N(age)=N(age-1)*sr(age-1);
end;
L = N(1:retage);
R = N(retage+1:maxage);

end


% OLG step in each iteration
% ------------------------------------------------------------------------
% given an interest rate "ret",
% Print implied assets and output and fval = rel - implied interest rate
% convergence => fval < eps
function [fval,ass,Y] = func_olg(ret)

global retage maxage L R N_age replrate tau sr lambda delta

% MPK = MC_k = interest rate + discount (from firm f.o.c.)
mpk = ret + delta;

% compute wage given MPK
wage = func_firm(mpk);

% compute pension contribution
tau = func_pens(L,R,replrate); % == 0 w/ current calibration

% compute income for each cohort
inc_det = func_inc(wage,tau,replrate,retage,maxage);

% compute gross saving from hh optimization problem
% (by age, shocks, and cash at hand) 
[sav_age, vfun] = func_hh_stoc(det_inc, ret, maxage, retage, sr, lambda);

% we want a unique saving level for each cohort:
    % @TODO determine life-cycle profile
    % @TODO aggregation by shock

% aggregation saving(cohort=t) = sav_t  * size_t
%ass = func_aggr(sav_age,N_age);

% compute capital level == saving => output (Y) implied by ret
[mpk,Y] = func_mpk(ass, L);

% update interest rate
retnew = mpk - delta; 

% change in interest rate
fval = ret - retnew;

end

% compute wages given mpk
% ------------------------------------------------------------------------
function wage=func_firm(mpk)

global alpha delta

k = (alpha/mpk)^(1/(1-alpha));
wage = (1-alpha)*k^alpha;

end


% pension system conribution/taxation? PAYGO?
% ------------------------------------------------------------------------
function tau = func_pens(L,R,replrate)

tau = replrate*R ./ (L + replrate * R);

end


% compute (deterministic) income profile over time given wage and pension  
% ------------------------------------------------------------------------
function inc=func_inc(wage,tau,replrate,retage,maxage);

inc=zeros(maxage,1);
inc(1:retage)=wage*(1-tau);
inc(retage+1:maxage)=replrate*wage*(1-tau);

end


% household optimization problem (deterministic)
% ------------------------------------------------------------------------
function gr_sav=func_hh(inc,ret,maxage,sr,lambda)

global betta tetta 

mpw = ones(maxage, 1);       % marginal propensity to consume out of total wealth
hk = zeros(maxage, 1);                 % human capital
RF = ones(maxage, 1)+ret;
annf = ones(maxage, 1);        % annuitization factor
for age = maxage-1 : -1 : 1,
    annf(age+1) = (sr(age) + lambda * (1-sr(age)))/sr(age);
    RF(age+1) = (1+ret) * annf(age+1);
    b = (betta * sr(age) * RF(age+1)^(1-tetta))^(-1/tetta);
    mpw(age) = [b * mpw(age+1)/(1 + b * mpw(age+1))]';

    % human capital: future discounted labor income
    hk(age) = (hk(age+1)+inc(age+1))/RF(age+1);
end;

% forward loop: optimal decisions:
gr_sav = zeros(maxage, 1);    % gr_sav = x-cons
net_sav = zeros(maxage, 1);   % net_sav = y-cons
x = zeros(maxage, 1);         % x = ass+y
cons = zeros(maxage, 1);      % consumption
tot_w = zeros(maxage, 1);     % total wealth (tot_w = ass+y+hk)
ass = zeros(maxage, 1);       % finanical assets
mpx = zeros(maxage, 1);

tot_w(1) = ass(1)+inc(1)+hk(1); 
for age = 1 : maxage,
    cons(age) = mpw(age) .* tot_w(age);
    x(age) = tot_w(age) - hk(age); 
    mpx(age)=cons(age)/x(age);
    ass(age) = x(age) - inc(age);
    gr_sav(age) = x(age) - cons(age);
    net_sav(age) = inc(age) - cons(age);
    if age < maxage,
        tot_w(age+1) = (tot_w(age)-cons(age))*RF(age+1);
    end;
end;
consgr=cons(2:end)./cons(1:end-1);

% check bugdet constraint:
df=ones(maxage,1);
af=annf(1);     % adjustment factor accounting for annuities
for age=1:maxage-1,
    af=af*annf(age+1);
    df(age+1)=1/(1+ret)^age/af;
end;
bc=sum(cons.*df-inc.*df);
if (abs(bc)>sqrt(eps))
    warning('present value budget constraint is violated');
end;

end


% household optimization problem (Stochastic)
% ------------------------------------------------------------------------
function [gr_sav, value]=func_hh_stoc(det_inc,ret,maxage,retage,sr,lambda)

global betta ret maxage retage sr lambda tetta opt_det opt_ny nshock nx curv grdfac

% Compute Markov Chain parameters
[nshock, pini, gridy,pi]=stoch_inc_proc(opt_det, opt_ny);

% grids and decisions rules:
gridx = zeros(maxage,nshock,nx); % years of life * possible distinct numbers for permanent shock realizations * possible cash-at-hand-states (three dimensional!)
gridsav = zeros(nx,1); % savings grid which we need in the alternative specification implemented here (one dimensional!)
gridass = zeros(maxage,nshock,nx); % to store associated assets => years of life * possible distinct numbers for permanent shock realizations * possible cash-at-hand-states (three dimensional!)
cfun = zeros(maxage,nshock,nx); % policy function in all states and ages for different cash-at-hand levels
vfun = zeros(maxage,nshock,nx); % value function in all states and ages for different cash-at-hand levels
vpfun = zeros(nx,1); % first derivative of the value function (w.r.t.consumption) (? why is the formula 8.62 in the dynmaic macro notes different?)
vptrans = zeros(maxage,nshock,nx); % transformed value function in all states and ages for different cash-at-hand levels

% savings grid: hold it constant:
maxsav=grdfac;
gridsav(2:nx)=makegrid(0.0,grdfac,nx-1,curv);
gridsav(1)=0.0;

% solve final period problem
% -----------------------------------------------------------------------
for yc=1:nshock
    
    % income at maxage (given deterministic income component)
    inc = det_inc(maxage)*gridy(yc);
    
    % set min grid
    minx=max(inc,sqrt(eps));
    
    % set max grid
    maxx=gridsav(nx)*(1.0+ret)+inc;
    
    % construct grid (given age and shock)
    gridx(maxage,yc,:)=linspace(minx,maxx,nx);
    
    % consumption = CAH
    cfun(maxage,yc,:)=gridx(maxage,yc,:);
    
    % assets from prevous period as implied dy CAH
    gridass(maxage,yc,:)=(gridx(maxage,yc,:)-inc)/(1+r);
    
    % value function
    vfun(maxage,yc,:)=U(cfun(maxage,yc,:));
    
    % marginal utility for Euler equation
    vpfun(:)=MUc(cfun(maxage,yc,:));
    vptrans(maxage,yc,:)=vpfun.^(-1.0/tetta);
end

% Age loop
% ------------------------------------------------------------------------
for jc = maxage-1:-1:1

    % Shock realization (today) loop
    % --------------------------------------------------------------------
    for yc=1:nshock
        
        % Cash at hand grid loop
        % ----------------------------------------------------------------
        for xc=2:nx
            
            vp = zeros(ny,1);
            
            % Solve tomorrows problem
            % ------------------------------------------------------------
           
            % loop over possible shock tomorrow
            for ycc=1:ny
                
                % income tomorrow (ASSUME == TODAY?)
                incp1 = det_inc(jc+1)*gridy(ycc);
                
                % Maximum cash on hand tomorrow:
                cah = max(sqrt(eps),incp1+(1.0+ret)*gridsav(xc));
                
                % Interpolate derivative of value function
                % --------------------------------------------------------
                if ( cah<gridx(jc+1,ycc,1))
                    disp('how can this be?')
                end
                
                if ( cah>gridx(jc+1,ycc,nx) )
                    % if out of bounds simply set it to decision at nx:
                    vptr = vptrans(jc+1,ycc,nx);
                else
                    vptr = interp1(squeeze(gridx(jc+1,ycc,:)),squeeze(vptrans(jc+1,ycc,:)),cah);
                end
                
                % compute value
                vp(ycc)=vptr.^(-tetta);
                
            end
            
            % Solve Euler equation
            expvp = betta*sr(jc)*(1.0+ret)*sum(pi(yc,:)*vp(:));
            cfun(jc,yc,xc)=invut(expvp);
            
            % update endogenous x-grid for cah:
            gridx(jc,yc,xc)=gridsav(xc)+cfun(jc,yc,xc);
        end
        
        % income in current period/age:
        inc = det_inc(jc)*gridy(yc);
        
        % decision at minx
        %-----------------------------------------------------------------
        minx=max(inc,sqrt(eps));
        
        % adjust x grid
        if (minx<gridx(jc,yc,2))
            gridx(jc,yc,1)=minx;
        else   
            gridx(jc,yc,1)=0.9*gridx(jc,yc,2);
        end
        
        % Compute optimal consumption for minx
        cfun(jc,yc,1)=gridx(jc,yc,1);
        
        % decision for each grid point
        %-----------------------------------------------------------------
        
        % assets implied by all x: 
        gridass(jc,yc,:)=(gridx(jc,yc,:)-inc)/(1+r); % <= TARGET
        
        % Update vfun and vpfun
        vpfun(:)=MUc(cfun(jc,yc,:));
        vptrans(jc,yc,:)=vpfun(:).^(-1.0/tetta);
        
        % compute value function
        for xc=1:nx
            
            v=zeros(nshock, 1);
            
            for ycc=1:ny
                
                % income tomorrow:
                incp1 = det_inc(jc+1)*gridy(ycc);
                
                % cah tomorrow
                cah = max(sqrt(eps),incp1+(1.0+ret)*gridsav(xc));
                
                % this should never be the case:
                if ((cah+0.0001)<gridx(jc+1,ycc,1)),
                    warning('How can this be ?');
                end;
                % linear interpolation:
                v(ycc)=func_intp(squeeze(gridx(jc+1,ycc,:)),squeeze(vfun(jc+1,ycc,:)),cah);
            end
            
            % update value function
            expv=sum(pi(yc,:)*v(:));
            
            % update Bellman equation
            vfun(jc,yc,xc)=U(cfun(jc,yc,xc))+betta*sr(jc)*expv;
            
            % print outputs
            % ------------------------------------------------------------
            
            % Saving grid == assets
            gr_sav = gridass;
            
            % value function 
            value = vfun;
    end
    
end
end


% Aggregation function
% ------------------------------------------------------------------------
function ass = func_aggr(sav_age,N_age)

ass = sum(sav_age.*N_age);

end


% MPK function
% ------------------------------------------------------------------------
function [mpk,Y] = func_mpk(ass, L)

global alpha

Y = ass.^alpha * L.^(1-alpha);
ky = ass./Y;
mpk = alpha * ky.^(-1);

end


% NEW - Markov chain income process
% ------------------------------------------------------------------------
function [ny, pini, gridy,pi]=stoch_inc_proc(opt_det, opt_ny)

if (opt_det==1),
    ny = 1;
    pini = 1.0;
    gridy = 1.0;
    pi = 1.0;
else
    
    if (opt_ny==1)
        % number of income shocks
        ny = 5;
        % transition probability
        rhoeta=0.98;
        % variance of "permanent" shock
        % taken from Campbell, Viceira, Ch. 7
        vareta=0.01;
        
        % Markov chain:
        % pi is the transition matrix of the Markov chain
        % gridy is the discretized state space of y_t
        [pi,gridy] = markovappr(rhoeta,sqrt(vareta),2,ny);

        
        % ii) exact by eigenvector decomposition:
        [v,d] = eig(pi');
        v = v./sum(v);
        pini = v(:,1);
        
        % take exponent and rescale income shocks such that mean is one:
        gridy=exp(gridy)';
        gridy = gridy / sum(gridy.*pini);
        

        
    else
        
        % Alternative -- taken from Kr?ger and Ludwig (2007) using only two
        % states
        ny = 2;
        
        % transition probability and variance
        rhoeta=0.97;
        vary=0.08;    % taken from Storesletten, Telmer, Yaron
        
        % shock
        epsil=sqrt(vary/(4.0*rhoeta*(1.0-rhoeta)));
        
        % Markov chain
        [pini,pi,gridy]=mchain(rhoeta,epsil);
    end
end
end
function [pini,pi,gridy]=mchain(rhoeta,epsil)

% Transition Probabilities
pi=rhoeta*ones(2,2);
pi(1,2)=1.0-rhoeta;
pi(2,1)=1.0-rhoeta;

% Initial Distribution
pini=0.5*ones(2,1);

gridy=zeros(2,1);
gridy(1)=exp(1.0-epsil);
gridy(2)=exp(1.0+epsil);
gridy=2.0*gridy/(sum(gridy));

end 


% NEW - Utility functions
% ------------------------------------------------------------------------
function u = U(c)
global tetta

% utility
if (abs(tetta-1-0)<sqrt(eps)),
    u = log(c);
else
    u = c.^(1.0-tetta)/(1.0-tetta);
end;
end     
function muc=MUc(c)
global tetta

% maringal utility
muc = c.^(-tetta);
end    
function invut=invut(marg)
global tetta

% invert utility for c
invut=marg.^(-1.0/tetta);
end    


