function lcpfchoice
% life-cycle portfolio choice model
% - solution of policy functions
% - MC simulation

close all
rho = 0.02;                                                                                  
tetta = 150.0;
pssi= 0.5
rf = 0.02;          % risk-free return
mu = 0.05;          % mean risky return
sig = 0.15;         % standard deviation of stock returns
ns = 10000;          % number of stochastic simulations
coh0 = 1.0;
inc0 = 1.0;         % net labor income during working period
rplr = 0.6;         % net replacement rate at retirement

nj = 61;            % maximum age (real age of 80)
jr = 46;            % retirement age (real age of 65) 

% transformation:
betta = 1/(1+rho);                                                          % define beta as 1/(1+rho)
inc = zeros(nj,1);                                                          % create zero matrix of size 61x1 for labor income??
inc(1:jr)=inc0;                                                             % assigning the first 46 values of vector inc0 to 1, where jr=46
inc(jr+1:nj)=rplr*inc0;                                                     % assigning from 47th element onwards to 61st the value of rplr*inc0(0.6*61)

% policy functions:
[mpc,hatalph]=func_pol(nj,betta,tetta,rf,mu,sig,pssi);                           % create policy function with inputs are: maximum age, preference discount factor, degree of risk aversion, risk-free return,mean risky return, standard deviation of stock returns.outputs: marginal propensity to consume, alphhat

% human capital
hk = zeros(nj,1);                                                            %create zero matrix of dimension 61*1 for human capital
for jc=nj-1:-1:1;
    hk(jc)=(hk(jc+1)+inc(jc+1))/(1+rf);                                     %calculate human capital for every successive periods
end;

% plot policy functions:
agevec=[1:nj]'+20-1;
figure;
pl=plot(agevec,mpc,'g-');
set(pl,'Linewidth',3);
title('marginal propensity to consume');

figure;
pl=plot(agevec,hatalph*ones(nj,1),'g-');
set(pl,'Linewidth',3);
title('transformed portfolio share');

% simulation
randn('state',0);        % initiazes the state of the random number generator
mret = zeros(nj,1);      % create zero matrix of size 61x1 for
mhatpfret = zeros(nj,1); % create zero matrix of size 61x1 for
mpfret = zeros(nj,1);    % create zero matrix of size 61x1 for
mcons = zeros(nj,1);     %  create zero matrix of size 61x1 for
malph = zeros(nj,1);     %  create zero matrix of size 61x1 for
mass = zeros(nj,1);      %  create zero matrix of size 61x1 for 
mtotwealth = zeros(nj,1);   %  create zero matrix of size 61x1 for mean total wealth
mcoh = zeros(nj,1);      %  create zero matrix of size 61x1 for
mtotsav = zeros(nj,1);    % create zero matrix of size 61x1 for mean total saving
mfinsav = zeros(nj,1);      %  create zero matrix of size 61x1 for mean financial saving
Rlong=[];
for sc=1:ns,
    lnR = randn(nj,1);      % log gross return (standardized)
    lnR = lnR*sig+log(1.0+mu);       % log gross return
    R = exp(lnR-0.5*sig^2);           % gross return
    ret = R-1.0;            % simple retur
    Rlong=[Rlong;R];
    
    [cons,consgr,coh,ass,totwealth,totsav,finsav,alph,hatpfret,pfret]=func_fwdsolve(coh0,hk,inc,nj,mpc,hatalph,rf,ret);
    
    if (sc==1),
        figure;
        pl=plot(agevec,ret,'g-');
        set(pl,'Linewidth',3);
        title('rate of return in first simulation');
        
        figure;
        pl=plot(agevec,cons,'g-');
        set(pl,'Linewidth',3);
        title('consumption in first simulation');

        figure;
        pl=plot(agevec,ass,'g-');
        set(pl,'Linewidth',3);
        title('assets in first simulation');

        figure;
        I = (alph(1:nj-1)>0.0);
        pl=plot(agevec(I),alph(I),'g-');
        set(pl,'Linewidth',3);
        title('share invested in risky asset in first simulation');

        figure;
        pl=plot(agevec(2:nj),pfret(2:nj),'g-');
        set(pl,'Linewidth',3);
        title('portfolio return in first simulation');
    end;
    
    % means:
    mret = mret + 1.0/ns*ret;
    mpfret = mpfret + 1.0/ns*pfret;
    mcons = mcons + 1.0/ns*cons;
    mass = mass + 1.0/ns*ass;
    mtotwealth = mtotwealth + 1.0/ns*totwealth;
    mcoh = mcoh + 1.0/ns*coh;
    malph(I) = malph(I) + 1.0/ns*alph(I);
    mtotsav = mtotsav + 1.0/ns*totsav;
    mfinsav = mfinsav + 1.0/ns*finsav;
end;

disp(['mean of R: ', num2str(mean(Rlong))])


 figure;
 pl=plot(agevec,mcons,'g-');
 set(pl,'Linewidth',3);
 title('mean consumption');

 figure;
 pl=plot(agevec,mass,'g-');
 set(pl,'Linewidth',3);
 title('mean assets');

 figure;
 pl=plot(agevec,mcoh,'g-');
 set(pl,'Linewidth',3);
 title('mean cash on hand');

 figure;
 pl=plot(agevec(1:nj-1),mcons(2:nj)./mcons(1:nj-1)-1.0,'g-');
 set(pl,'Linewidth',3);
 title('mean consumption growth');

 figure;
 I = (malph(1:nj-1)>0.0);
 pl=plot(agevec(I),malph(I),'g-');
 set(pl,'Linewidth',3);
 title('mean share invested in risky asset');

 figure;
 pl=plot(agevec(2:nj),mpfret(2:nj),'g-');
 set(pl,'Linewidth',3);
 title('mean portfolio return');

 disp('THE END');
 beep; beep; beep;
 % -------------------------------------------------------------------

end
% -------------------------------------------------------------------

% -------------------------------------
function [mpc,hatalph]=func_pol(nj,betta,tetta,rf,mu,sig,pssi)

lnmu=log(1.0+mu);
mpc = ones(nj,1);                  % marginal propensity to consume out of cash on hand
hatalph=(lnmu-log(1.0+rf)+sig^2/2)/(tetta*sig^2);   % approximate solution to portfolio share
mup = hatalph*lnmu + (1-hatalph)*log(1+rf) + 0.5*hatalph*(1.0-hatalph)*sig^2;
sigp = hatalph^2*sig^2;
for jc = nj-1 : -1 : 1,
    b = mpc(jc+1)^(-1) * (betta * exp( (1.0-1/pssi)*(mup+(1.0-tetta)*sigp^2/2) ) )^(pssi);  
    mpc(jc) = 1.0/(1.0 + b);
end;

end
% ---------------------------------


% ---------------------------------
function [cons,consgr,coh,ass,totwealth,totsav,finsav,alph,hatpfret,pfret]=func_fwdsolve(coh0,hk,inc,nj,mpc,hatalph,rf,ret)

totsav = zeros(nj, 1);       % savings = coh-c
finsav = zeros(nj, 1);       % savings = coh-c
coh = zeros(nj, 1);         % coh = ass*R
cons = zeros(nj, 1);      % consumption
ass = zeros(nj, 1);       % finanical assets
totwealth = zeros(nj, 1);
alph = zeros(nj,1);
pfret =zeros(nj,1);
hatpfret =zeros(nj,1);

coh(1)=coh0;
totwealth(1) = coh0+hk(1);  % initial total wealth
for jc = 1 : nj,
    % consumption and cash on hand
    cons(jc)=mpc(jc)*totwealth(jc);
    coh(jc)=totwealth(jc)-hk(jc);
    if (jc==1),
        hatpfret(jc)=rf+hatalph*(ret(jc)-rf);
    end;
    ass(jc)=coh(jc)-inc(jc);
    
    % savings
    totsav(jc)=totwealth(jc)-cons(jc);
    finsav(jc)=coh(jc)-cons(jc);
    
    if jc < nj,
        alph(jc)=hatalph*totsav(jc)/finsav(jc);
        hatpfret(jc+1)=rf+hatalph*(ret(jc+1)-rf);
        pfret(jc+1)=rf+alph(jc)*(ret(jc+1)-rf);
        totwealth(jc+1) = totsav(jc)*(1+hatpfret(jc+1));
    end;
end;
consgr=cons(2:end)./cons(1:end-1);

end
% ---------------------------------






