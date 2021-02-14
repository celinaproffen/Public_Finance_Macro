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

maxage = 80;        % corresponds to real age 100
retage = 45;        % corresponds to actual age 65
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
% ------------------------------------------------------------------------

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
% ------------------------------------------------------------------------

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
inc = func_inc(wage,tau,replrate,retage,maxage);

% compute gross saving from hh (by age) optimization problem
sav_age = func_hh(inc,ret,maxage,sr,lambda);

% aggregation saving(cohort=t) = sav_t  * size_t
ass = func_aggr(sav_age,N_age);

% compute capital level == saving => output (Y) implied by ret
[mpk,Y] = func_mpk(ass, L);

% update interest rate
retnew = mpk - delta; 

% change in interest rate
fval = ret - retnew;

end
% ------------------------------------------------------------------------

% compute wages given mpk
% ------------------------------------------------------------------------
function wage=func_firm(mpk)

global alpha delta

k = (alpha/mpk)^(1/(1-alpha));
wage = (1-alpha)*k^alpha;

end
% ------------------------------------------------------------------------

% pension system conribution/taxation? PAYGO?
% ------------------------------------------------------------------------
function tau = func_pens(L,R,replrate)

tau = replrate*R ./ (L + replrate * R);

end
% ------------------------------------------------------------------------

% compute income profile over time given wage and pension  
% ------------------------------------------------------------------------
function inc=func_inc(wage,tau,replrate,retage,maxage);

inc=zeros(maxage,1);
inc(1:retage)=wage*(1-tau);
inc(retage+1:maxage)=replrate*wage*(1-tau);

end
% ------------------------------------------------------------------------

% household optimization problem
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
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
function ass = func_aggr(sav_age,N_age)

ass = sum(sav_age.*N_age);

end
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
function [mpk,Y] = func_mpk(ass, L)

global alpha

Y = ass.^alpha * L.^(1-alpha);
ky = ass./Y;
mpk = alpha * ky.^(-1);

end
% ------------------------------------------------------------------------








