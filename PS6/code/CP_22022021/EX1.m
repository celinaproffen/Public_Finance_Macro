% exercize
function exercise

close all

global alpha betta tetta delta replrate tau lambda 
global retage maxage L R N N_age sr nx curv grdfac 

rho = 0.05;  % 0.01
delta = 0.05;
alpha = 0.33; % 0.33; % 0.1;
ret = 0.07;
% tetta = 1.0;
tetta=1.6, % this is our new tetta for comparison in problem 2
betta = 1/(1+rho);
lambda = 1.0;
opt_PAYG=0
opt_ra=0
maxage = 80;        % corresponds to real age 100 == nj
retage = 45;        % corresponds to actual age 65 == jr
%replrate = 0.0;     % ?????
replrate = 0.6;     % This is the new rr=0.6 as suggested in problem 2
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
opt_det = 1;       
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
for i = 1:n
    ret = retvec(i);
    [delta_ret, assvec(i)] = func_olg(ret);
end

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

% iterate until r* s.t. => delta_ret < tol 
for it=1:maxit
    
    % given ret, compute differene between demand and supply r
    [delta_ret, ass, Y] = func_olg(ret);
    
    % if delta_ret(ret) < tol, STOP 
    if (abs(delta_ret)<tol)
        disp('convergence');
        break;
        
    % delta_ret(ret) > tol
    else
        
        % update ret using df as a smoothing parameter
        ret = ret - df * delta_ret;
        disp(['iteration #', num2str(it)]);
        disp(['guess of rate of return: ', num2str(ret)]);
    end
end

if (it>=maxit)
    warning('no convergence');
end
disp(['equilibrium interest rate: ', num2str(ret)]);
disp(['equilibrium contribution rate: ', num2str(tau)]);
disp(['equilibrium capital output ratio: ', num2str(ass/Y)]);

end


% compute the size of each cohort given min, max and survival rate sr
% ------------------------------------------------------------------------
function [N,L,R]=func_pop(sr,maxage,retage)

N = zeros(maxage,1);
L = zeros(maxage,1);
R = zeros(maxage,1);
N(1) = 100;
for age=2:maxage
    N(age)=N(age-1)*sr(age-1);
end
L = N(1:retage);
R = N(retage+1:maxage);

end


% OLG step in each iteration
% ------------------------------------------------------------------------
% given an interest rate "ret",
% Print implied assets and output and delta_ret = rel - implied interest rate
% convergence => delta_ret < eps
function [delta_ret,ass,Y] = func_olg(ret)

global retage maxage L R replrate tau sr  delta gridy pini pi

% MPK = MC_k = interest rate + discount (from firm f.o.c.)
mpk = ret + delta;

% compute wage given MPK
wage = func_firm(mpk);

% compute pension contribution
%NEW! equivalent income variation analysis
tau = func_pens(L,R,replrate); % == 0 w/ current calibration
for opt_PAYG=0:1
    if opt_PAYG==0
        replrate=0.0
        tau = func_pens(L,R,replrate);
        det_inc = func_inc(wage,tau,replrate,retage,maxage);
        [gridx,gridsav,gridass,cfun,vfun]= func_hh_stoc(det_inc,ret);
        V=vfun(0,:,:)
    else
        replrate=0.6
        tau = func_pens(L,R,replrate); 
        det_inc = func_inc(wage,tau,replrate,retage,maxage);
        [gridx,gridsav,gridass,cfun,vfun]= func_hh_stoc(det_inc,ret)
        V_bar=vfun(0,:,:)
    end
end
%Equivalent income variation   
g=eqiv(V,V_bar,tetta)

%NEW- variation across different values of tetta 
for opt_ra=0:1
    if opt_ra==0
       tetta=1.6
       g=eqiv(V,V_bar,tetta)
    else 
       tetta=1.2
       g=eqiv(V,V_bar,tetta)
    end
end
% compute deterministic component of income for each cohort
det_inc = func_inc(wage,tau,replrate,retage,maxage);

% compute gross saving from hh optimization problem
% (by age, shocks, and cash at hand) 
[gridx,gridsav,gridass,cfun,vfun] = func_hh_stoc(det_inc,ret);

% aggregation 
[Phi,PhiAss,ass] = func_aggr_stoc(gridx, gridsav, cfun, det_inc);

% compute capital level == saving => output (Y) implied by ret
[mpk,Y] = func_mpk(ass, L);

% update interest rate
retnew = mpk - delta; 

% change in interest rate
delta_ret = ret - retnew;

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
function inc = func_inc(wage,tau,replrate,retage,maxage)

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
for age = maxage-1 : -1 : 1
    annf(age+1) = (sr(age) + lambda * (1-sr(age)))/sr(age);
    RF(age+1) = (1+ret) * annf(age+1);
    b = (betta * sr(age) * RF(age+1)^(1-tetta))^(-1/tetta);
    mpw(age) = [b * mpw(age+1)/(1 + b * mpw(age+1))]';

    % human capital: future discounted labor income
    hk(age) = (hk(age+1)+inc(age+1))/RF(age+1);
end

% forward loop: optimal decisions:
gr_sav = zeros(maxage, 1);    % gr_sav = x-cons
net_sav = zeros(maxage, 1);   % net_sav = y-cons
x = zeros(maxage, 1);         % x = ass+y
cons = zeros(maxage, 1);      % consumption
tot_w = zeros(maxage, 1);     % total wealth (tot_w = ass+y+hk)
ass = zeros(maxage, 1);       % finanical assets
mpx = zeros(maxage, 1);

tot_w(1) = ass(1)+inc(1)+hk(1); 
for age = 1 : maxage
    cons(age) = mpw(age) .* tot_w(age);
    x(age) = tot_w(age) - hk(age); 
    mpx(age)=cons(age)/x(age);
    ass(age) = x(age) - inc(age);
    gr_sav(age) = x(age) - cons(age);
    net_sav(age) = inc(age) - cons(age);
    if age < maxage
        tot_w(age+1) = (tot_w(age)-cons(age))*RF(age+1);
    end
end
consgr=cons(2:end)./cons(1:end-1);

% check bugdet constraint:
df=ones(maxage,1);
af=annf(1);     % adjustment factor accounting for annuities
for age=1:maxage-1
    af=af*annf(age+1);
    df(age+1)=1/(1+ret)^age/af;
end
bc=sum(cons.*df-inc.*df);
if (abs(bc)>sqrt(eps))
    warning('present value budget constraint is violated');
end

end


% NEW! household optimization problem (stochastic)
% ------------------------------------------------------------------------
function [gridx,gridsav,gridass,cfun,vfun]=func_hh_stoc(det_inc,ret)

global betta tetta opt_det opt_ny nshock nx curv grdfac pi gridy pini maxage sr

% Compute Markov Chain parameters
[nshock, pini, gridy, pi]=stoch_inc_proc(opt_det, opt_ny);

% grids and decisions rules:
gridx = zeros(maxage,nshock,nx); % years of life * possible distinct numbers for permanent shock realizations * possible cash-at-hand-states (three dimensional!)
gridsav = zeros(nx,1); % savings grid which we need in the alternative specification implemented here (one dimensional!)
gridass = zeros(maxage,nshock,nx); % to store associated assets => years of life * possible distinct numbers for permanent shock realizations * possible cash-at-hand-states (three dimensional!)
cfun = zeros(maxage,nshock,nx); % policy function in all states and ages for different cash-at-hand levels
vfun = zeros(maxage,nshock,nx); % value function in all states and ages for different cash-at-hand levels
vpfun = zeros(nx,1); % first derivative of the value function (w.r.t.consumption) (? why is the formula 8.62 in the dynmaic macro notes different?)
vptrans = zeros(maxage,nshock,nx); % transformed value function in all states and ages for different cash-at-hand levels

% savings grid: hold it constant:
maxsav = grdfac;
gridsav(2:nx)=makegrid(0.0, maxsav, nx-1, curv);
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
    gridass(maxage,yc,:)=(gridx(maxage,yc,:)-inc)/(1+ret);
    
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
            
            vp = zeros(nshock, 1);
            
            % Solve tomorrows problem
            % ------------------------------------------------------------
           
            % loop over possible shock tomorrow
            for ycc=1:nshock
                
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
        gridass(jc,yc,:)=(gridx(jc,yc,:)-inc)/(1+ret);    % <= TARGET
        
        % Update vfun and vpfun
        vpfun(:)=MUc(cfun(jc,yc,:));                    % <= TARGET
        vptrans(jc,yc,:)=vpfun(:).^(-1.0/tetta);
        
        % compute value function
        for xc=1:nx
            
            v=zeros(nshock, 1);
            
            for ycc=1:nshock
                
                % income tomorrow:
                incp1 = det_inc(jc+1)*gridy(ycc);
                
                % cah tomorrow
                cah = max(sqrt(eps),incp1+(1.0+ret)*gridsav(xc));
                
                % this should never be the case:
                if ((cah+0.0001)<gridx(jc+1,ycc,1))
                    warning('How can this be ?');
                end
                % linear interpolation:
                v(ycc)=func_intp(squeeze(gridx(jc+1,ycc,:)),squeeze(vfun(jc+1,ycc,:)),cah);
            end
            
            % update value function
            expv=sum(pi(yc,:)*v(:));
            
            % update Bellman equation
            vfun(jc,yc,xc)=U(cfun(jc,yc,xc))+betta*sr(jc)*expv;
            
        end
    
    end
end
% Auxiliary functions
% ------------------------------------------------------------------------

% Interpolation 
    function fv = func_intp(x,func,xp)
        
        
        n = length(x);
        if ( xp>x(n) ),
            % fv = func(n);
            fv=func_extrapol(x(n-1),x(n),func(n-1),func(n),xp);
        elseif (xp<x(1)),
            % fv = func(1);
            fv=func_extrapol(x(1),x(2),func(1),func(2),xp);
        else
            fv = interp1(x,func,xp);
        end;
        
    end
% extrapolation function
    function y=func_extrapol(x1,x2,y1,y2,x)
        
        % simple linear extrapolation
        
        m = (y2-y1)/(x2-x1);
        y = y1 + m*(x-x1);
        
    end

% makegrid
    function grd = makegrid(x1,x2,n,c)
    % makes curved grid according to curvature parameter c
    scale=x2-x1;
    grd(1)=x1;
    grd(n)=x2;
    for i=2:n-1
        grd(i)=x1+scale*((i-1.0)/(n-1.0))^c;
    end
end 
end 


% Aggregation function (deterministic)
% ------------------------------------------------------------------------
function ass = func_aggr(sav_age,N_age)

ass = sum(sav_age.*N_age);

end


% NEW! Aggregation function (stochastic)
% ------------------------------------------------------------------------
function [Phi,PhiAss,ass]=func_aggr_stoc(gridx, gridsav, cfun, det_inc)
    
global nshock ret maxage nx pi gridy sr pini N
    
disp('aggregation and cross-sectional measure');
    
% compute total population
totpop = sum(N);
frac = N./totpop;
    
% Compute Cross sectional distributions and aggregate variables
% --------------------------------------------------------------------
    
% distribution of assets conditional by age and shock
Phi = zeros(maxage, nshock, nx);
    
% distribution of assets
PhiAss = zeros(nx,1); 
    
% Distribution of newborns over cash at hand
for yc = 1:nshock
    
    % income in current period/age:
    inc = det_inc(1)*gridy(yc);
    
	% initial cash-on-hand:
    cahini = inc;
        
	% find values and position of positive 
	[vals,inds]=basefun(gridx(1,yc,:),cahini,nx);
        
	% compute transition 
    Phi(1,yc,inds(1))=vals(1)*pini(yc)*frac(1);
	Phi(1,yc,inds(2))=vals(2)*pini(yc)*frac(1);
 
end
    
% compute transfer function
%---------------------------------------------------------------------
for jc=2:maxage
        
	% transfer function
    TT = zeros(nshock, nx, nshock, nx);
        
    for xc=1:nx
    	for yc=1:nshock
        	for ycc=1:nshock
                
            	% income (wages and pensions) in current period/age:
                inc = det_inc(jc)*gridy(ycc);
                
                % cash on hand: x=a*(1+r)+y = s(-1)*(1+r)+y;
                cah = inc+ (1.0+ret)*gridsav(xc);
                
                [vals,inds]=basefun(gridx(jc,ycc,:),cah,nx);
                
                TT(ycc,inds(1),yc,xc)=vals(1)*pi(yc,ycc);
                TT(ycc,inds(2),yc,xc)=vals(2)*pi(yc,ycc);
            end  
        end    
    end
end
    
for xc=1:nx
    for yc=1:nshock
        for xcc=1:nx
            for ycc=1:nshock
            	% transfer distribution:
                Phi(jc,ycc,xcc)=Phi(jc,ycc,xcc)+Phi(jc-1,yc,xc)*TT(ycc,xcc,yc,xc)*sr(jc-1);
            end
        end
    end
end

% Check if distribution sums to 1
sumprob=sum(sum(sum(Phi(:,:,:))));
    if ( ( sumprob < 0.999 ) || ( sumprob > 1.001) )
        beep; beep; beep;
        warning('distribution fucked up');
    end
    
% Check if Grid is Big enough
sumprob=sum(sum(Phi(:,:,nx)));
if (sumprob > 0.001 )
    beep; beep; beep;
    warning('grid too small -- increase your grid');
    pause
end

% NOTE: KEEP RET?

% aggregation
%-------------------------------------------------------------------------
for jc=1:maxage
    for yc=1:nshock
        for xc=1:nx
            
            PhiAss(xc)=PhiAss(xc)+Phi(jc,yc,xc);
            
            % asset holdings = capital stock in general equilibrium
            ass = totpop*Phi(jc,yc,xc)*gridsav(xc);
            
            % consumption level
            cons = totpop*Phi(jc,yc,xc)*cfun(jc,yc,xc);
            
            %lab = totpop*Phi(jc,yc,xc)*gridy(yc)*epsi(jc);
            %ret=ret+totpop*Phi(jc,yc,xc)*gridy(yc)*(1.0-epsi(jc));
        end
    end
end

% Basefun
%-------------------------------------------------------------------------
 function [vals,inds]=basefun(grid_x,x,nx)
        % this subroutine returns the values and the indices of the two basis
        % functions that are positive on a given x in the grid_x
        
        % MF function to lookup the current position
        % I think this gives you the closest level in grid_x to x (here
        % cah) in term of the index of grid_x where to find that cah-level
        i = lookup(grid_x,x,0);
        
        disp([i, nx, "end"])
        
        if ((i+1) > nx)
            inds(1) = nx;
            inds(2) = nx;
            vals(2) = 0.0;
            vals(1) = 1.0;
        elseif (i == 0)
            inds(1) = 1;
            inds(2) = 1;
            vals(1) = 1.0;
            vals(2) = 0.0;
        else
            inds(1) = i;
            inds(2) = i+1;
            dist = grid_x(i+1)-grid_x(i);
            vals(2) = ( x-grid_x(i) )/dist;
            vals(1) = ( grid_x(i+1)-x )/dist;
        end
        
 end 

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
function [nshock, pini, gridy,pi]=stoch_inc_proc(opt_det, opt_ny)

if (opt_det==1)
    nshock = 1;
    pini = 1.0;
    gridy = 1.0;
    pi = 1.0;
else
    
    if (opt_ny==1)
        % number of income shocks
        nshock = 5;
        % transition probability
        rhoeta=0.98;
        % variance of "permanent" shock
        % taken from Campbell, Viceira, Ch. 7
        vareta=0.01;
        
        % Markov chain:
        % pi is the transition matrix of the Markov chain
        % gridy is the discretized state space of y_t
        [pi,gridy] = markovappr(rhoeta,sqrt(vareta),2,nshock);

        
        % ii) exact by eigenvector decomposition:
        [v,d] = eig(pi');
        v = v./sum(v);
        pini = v(:,1);
        
        % take exponent and rescale income shocks such that mean is one:
        gridy=exp(gridy)';
        gridy = gridy / sum(gridy.*pini);
        

        
    else
        
        % Alternative -- taken from Krueger and Ludwig (2007) using only two
        % states
        nshock = 2;
        
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
if (abs(tetta-1-0)<sqrt(eps))
    u = log(c);
else
    u = c.^(1.0-tetta)/(1.0-tetta);
end
end     
function muc=MUc(c)
global tetta

% marginal utility
muc = c.^(-tetta);
end    
function invut=invut(marg)
global tetta

% invert utility for c
invut=marg.^(-1.0/tetta);
end   
% NEW - Equivalent variation
function eqiv=EQIV(V,V_bar,tetta)
eqiv=(V/V_bar)^(1/(1-tetta))-1
end



