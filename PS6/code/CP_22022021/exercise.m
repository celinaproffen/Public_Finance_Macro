% exercise

% Comment celina: I did not entirely understand how this loops over
% returns... it uses the func_olg

function exercise

close all

global alpha betta tetta delta replrate tau lambda 
global retage maxage L R N_age sr 

rho = 0.05;  % 0.01
delta = 0.05;
alpha = 0.33; % 0.33; % 0.1;
ret = 0.09;
tetta = 1.0;
betta = 1/(1+rho);
lambda = 1.0; % ???

maxage = 80;    % corresponds to real age 100
retage = 45;    % corresponds to actual age 65
replrate = 0.0; % Retirement plan rate??

tol = 1e-4;
maxit = 100;
df = 0.1;

mr = readfile([],'MR.txt',3);
sr = 1.0-mr(21:end,1);

[N_age,L_age,R_age] = func_pop(sr,maxage,retage);
N = sum(N_age); % overall people alive ?
L = sum(L_age); % labor?
R = sum(R_age); % retirement?

% graph of A(r) and K(r)
retold = ret;

% graph of A(r)
retvec = [-delta+0.01:0.001:0.1]'; % vector of return grid points
n = length(retvec);
assvec = zeros(n,1); % vector of assets
for i = 1:n,
    ret = retvec(i);
    [fval, assvec(i)] = func_olg(ret);
end;

% graph of K(r):
kvec = (alpha./(retvec+delta)).^(1/(1-alpha)).*L;

figure; hold on;
plot(retvec,assvec,'b-','LineWidth',3);
plot(retvec,kvec,'r--','LineWidth',3);
legend('demand', 'supply');
title('supply and demand in OLG economy');
print -depsc equil.eps


% solution of OLG model:
ret = retold;

for it=1:maxit,
    [fval, ass, Y] = func_olg(ret);
    if (abs(fval)<tol),
        disp('convergence');
        break;
    else,
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
% ---------------------------------


% ---------------------------------
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
% ---------------------------------


% ----------------------------------------
function [fval,ass,Y] = func_olg(ret)

global retage maxage L R N_age replrate tau sr lambda delta

mpk = ret + delta; % marginal product of capital

wage = func_firm(mpk); % marginal product of labor that depends on capital

tau = func_pens(L,R,replrate); % ??? maybe the labor tax rate? depends on number of years at labor and retired

inc = func_inc(wage,tau,replrate,retage,maxage); % labor income

sav_age = func_hh(inc,ret,maxage,sr,lambda); % savings decisions giving rate of return (these are the gross savings)

ass = func_aggr(sav_age,N_age);

[mpk,Y] = func_mpk(ass, L); 

retnew = mpk - delta; 

fval = ret - retnew; % Change in rate of return - diFference in Value?

end
% ----------------------------------------


% -------------------------------------------
function wage=func_firm(mpk)

global alpha delta

k = (alpha/mpk)^(1/(1-alpha));
wage = (1-alpha)*k^alpha;

end
% -------------------------------------------

% ----------------------------------------------
function tau = func_pens(L,R,replrate)

tau = replrate*R ./ (L + replrate * R);

end
% ----------------------------------------------

% ---------------------------------------------
function inc=func_inc(wage,tau,replrate,retage,maxage);

% This is actually way more complicated.. We gave people a potentially
% risky (Markov process) income process, that had a deterministic age
% component and shock components (both idiosyncratic and aggregate, 
% although we lumped them together)

inc=zeros(maxage,1);
inc(1:retage)=wage*(1-tau);
inc(retage+1:maxage)=replrate*wage*(1-tau);

end
% ---------------------------------------------

% -------------------------------------
function gr_sav=func_hh(inc,ret,maxage,sr,lambda)

global betta tetta 

mpw = ones(maxage, 1);       % marginal propensity to consume out of total wealth
hk = zeros(maxage, 1);       % human capital
RF = ones(maxage, 1)+ret;    % Gross return
annf = ones(maxage, 1);        % annuitization factor
for age = maxage-1 : -1 : 1,
    annf(age+1) = (sr(age) + lambda * (1-sr(age)))/sr(age); % ?? 
    % NOTE: right now lambda = 1, hence annf = 1/how likely it was to become this old given you were age-1 in the last period
    RF(age+1) = (1+ret) * annf(age+1); % expected rate of return, accounting for the possibility of dying tomorrow
    b = (betta * sr(age) * RF(age+1)^(1-tetta))^(-1/tetta);
    mpw(age) = [b * mpw(age+1)/(1 + b * mpw(age+1))]'; % In the last period
    % you will consume evrything, before less (always accounting for 
    % possibility of dying, rate of return, time discount factor and RRA)

    % human capital: future discounted labor income
    hk(age) = (hk(age+1)+inc(age+1))/RF(age+1);
end;

% forward loop: optimal decisions:
gr_sav = zeros(maxage, 1);    % gr_sav = x-cons (wealth or cah - cons)
net_sav = zeros(maxage, 1);   % net_sav = y-cons (labor income - cons)
x = zeros(maxage, 1);         % x = ass+y
cons = zeros(maxage, 1);      % consumption
tot_w = zeros(maxage, 1);     % total wealth (tot_w = ass+y+hk)
ass = zeros(maxage, 1);       % financial assets
mpx = zeros(maxage, 1);       % MPC out of cah

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
% -------------------------------------

% ----------------------------------------
function ass = func_aggr(sav_age,N_age)

ass = sum(sav_age.*N_age);

end
% ----------------------------------------

% ----------------------------------------
function [mpk,Y] = func_mpk(ass, L)

global alpha

Y = ass.^alpha * L.^(1-alpha);
ky = ass./Y;
mpk = alpha * ky.^(-1);

end
% ----------------------------------------








