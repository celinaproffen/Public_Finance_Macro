% ++++++++++beg++++++++++++++++++++++++++++++++++++++++++++
function towards_olg

close all

global nj ny

tic

opt_det=false;          % 1=deterministic model
opt_nosr=0;         % 1=no survival risk
opt_ny = 2;             % 1=Markov chain with number of states, ny=5,
% 2=Markov chain with ny=2 (Krüger-Ludwig
% calibration)

% -------------------------------------------------------------------------
% SOLUTION

% calibration
func_calibr(opt_det,opt_nosr,opt_ny);

% solution of household model
[gridx,gridsav,gridass,cfun,vfun] = func_hh;

% aggregation
[Phi,PhiAss,ass] = func_aggr(gridx,gridsav,cfun,gridass);

% average life-cycle profiles
[labinclife,inclife,asslife,conslife,vallife] = lcprofile(Phi,gridass,cfun,vfun);
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% PLOTS

% plot of consumption policy for seleted ages and grid points and
% current shock state fix((ny+1)/2):
avec=[1,21,41,61];

figure;
pl=plot(squeeze(gridx(1:10))',squeeze(cfun(avec,fix((ny+1)/2),[1:10]))');
legend('age 20','age 40','age 60','age 80');
set(pl,'LineWidth',2);
title('consumption policy at different ages');
print ('-depsc', ['conspol_', num2str(opt_det), '_', num2str(opt_nosr), '.eps']);

% distribution by age and shock
for age = avec,
    figure;
    pl=plot(gridsav,squeeze(Phi(age,:,:))');
    set(pl,'LineWidth',4);
    title(['distribution (assets tomorrow) at age ', num2str(age+20-1)]);
    print ('-depsc', ['distr_age', num2str(age), '_', num2str(opt_det), '_', num2str(opt_nosr), '.eps']);
end;

% asset distribution
figure;
pl=plot(gridsav,PhiAss);
set(pl,'LineWidth',4);
title('asset distribution (assets tomorrow)');
print ('-depsc', ['PhiAss_', num2str(opt_det), '_', num2str(opt_nosr), '.eps']);

% plots of average life-cycle profiles
age=[20:nj+20-1];
figure;
pl=plot(age,labinclife); title('labor income'); xlabel('age');
set(pl,'LineWidth',2);
print ('-depsc', ['labinc_', num2str(opt_det), '_', num2str(opt_nosr), '.eps']);

figure;
pl=plot(age,inclife);  title('income'); xlabel('age');
set(pl,'LineWidth',2);
print ('-depsc', ['inc_', num2str(opt_det), '_', num2str(opt_nosr), '.eps']);

figure;
pl=plot(age,asslife); title('assets'); xlabel('age');
set(pl,'LineWidth',2);
print ('-depsc', ['ass_', num2str(opt_det), '_', num2str(opt_nosr), '.eps']);

figure;
pl=plot(age,conslife); title('consumption'); xlabel('age');
set(pl,'LineWidth',2);
print ('-depsc', ['cons_', num2str(opt_det), '_', num2str(opt_nosr), '.eps']);

cgr = conslife(2:end)./conslife(1:end-1);
figure;
pl=plot(age(1:end-1),cgr); title('consumption growth'); xlabel('age');
set(pl,'LineWidth',2);
print ('-depsc', ['consgr_', num2str(opt_det), '_', num2str(opt_nosr), '.eps']);

figure;
pl=plot(age,inclife-conslife); title('savings'); xlabel('age');
set(pl,'LineWidth',2);
print ('-depsc', ['sav_', num2str(opt_det), '_', num2str(opt_nosr), '.eps']);

figure;
pl=plot(age,vallife); title('value'); xlabel('age');
set(pl,'LineWidth',2);
print ('-depsc', ['value_', num2str(opt_det), '_', num2str(opt_nosr), '.eps']);
% -------------------------------------------------------------------------

disp(['time elapsed: ', num2str(toc)]);

end     % end function func_main
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function func_calibr(opt_det,opt_nosr,opt_ny)

global betta tetta r nj jr nx ny pi gridy netw pens sr epsi curv pini frac pop totpop grdfac

close all

r = 0.04;
rho = 0.04;
betta = 1/(1+rho);
tetta = 2;

nj=80;
jr=45;

nx=30;         % # of grid-points
curv=3.0;       % curvature of grid
grdfac=40;      % scaling factor of saving grid

% deterministic income component:
netw=1.0;
pens=0.4;
epsi=ones(nj,1);
if (jr<nj),
    epsi(jr+1:nj)=0.0;
end;

% survival rates
if opt_nosr,
    sr = ones(nj,1);
else
    mr = readfile([],'MR.txt',3);
    sr = 1.0-mr(21:21+nj-1,1);
end;

% population and fraction living in year...
pop=zeros(nj,1);
pop(1)=100;
for jc=2:nj,
    pop(jc)=pop(jc-1)*sr(jc-1);
end;
totpop=sum(pop);

% normalize population to one:
pop=pop/totpop;
totpop=1.0;
frac=pop./totpop;

% # of income states
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
        [pi,gridy] = markovappr(rhoeta,sqrt(vareta),2,ny);
        
        % compute invariant distribution
        % i) by simulation
%         pini = 1/ny*ones(ny,1);
%         for tc=1:10000,
%             pini = pi'*pini;
%         end;
        
        % ii) exact by eigenvector decomposition:
        [v,d] = eig(pi');
        v = v./sum(v);
        pini = v(:,1);
        
        % take exponent and rescale income shocks such that mean is one:
        gridy=exp(gridy)';
        gridy = gridy / sum(gridy.*pini);
        
    else
        
        % Alternative -- taken from Krüger and Ludwig (2007) using only two
        % states
        ny = 2;
        
        % transition probability and variance
        rhoeta=0.97;
        vary=0.08;    % taken from Storesletten, Telmer, Yaron
        
        % shock
        epsil=sqrt(vary/(4.0*rhoeta*(1.0-rhoeta)));
        
        % Markov chain
        [pini,pi,gridy]=mchain(rhoeta,epsil);
    end;
    
end;


end     % end function func_calibr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [gridx,gridsav,gridass,cfun,vfun] = func_hh

global betta tetta r jr nj nx ny pi gridy netw pens sr epsi grdfac curv

disp('solution of household model');

% grids and decisions rules:
gridx = zeros(nx,1);
gridass = zeros(nj,ny,nx);
cfun = zeros(nj,ny,nx);
vfun = zeros(nj,ny,nx);
vpfun = zeros(nj,ny,nx);

% (2) define cash-on-hand grid
%-------------------------------------------------------------------

% since borrowing costraint is a >= 0, the minimum cash-on-hand is
minx = sqrt(eps);

% the maximun is given by maximum possible income over the lifetime:
% maxx = max (labor) income + max (pension) income 
maxx = grdfac;

% construct grid with equidistant points
gridx(1:nx) = linspace(minx,maxx,nx);

% saving grid?
gridsav(1:nx)= linspace(0.0,maxx,nx);

% (3) solve problem at final period (nj)
%-------------------------------------------------------------------

for yc=1:ny     % for each possible shock in nj
    
    % Final period income 
    inc = epsi(nj)*netw*gridy(yc) + (1-epsi(nj))*pens;
    
    % Final period Consumption function = cash at hand (for all shocks)
    cfun(nj,yc,:) = gridx;
    
    % asset holdings, value function, including derivative
    gridass(nj,yc,:) = (gridx-inc)/(1+r);
    vfun(nj,yc,:) = U(cfun(nj,yc,:));
    vpfun(nj,yc,:) = MUc(cfun(nj,yc,:));
    
end

% Iterate Backwards
%-------------------------------------------------------------------

% iterate for each period/age
for jc=nj-1:-1:1
    
    % iterate for each possible shock today
    for yc=1:ny
        
        % iterate for each point cash at hand
        for xc=1:nx
            
            % compute expected value of future value function V_x'(x')
            vp = zeros(2,1);
            
            % for each possible shock tomorrow
            for ycc=1:ny
                
                % income tomorrow:
                incp1 = epsi(jc+1)*netw*gridy(ycc)+(1-epsi(jc+1))*pens;
                
                % compute x' given x today (c ??????) 
                cah = max(sqrt(eps),incp1+(1.0+r)*gridass(jc+1,ycc,xc));
                
                % Interpolate derivative of value function
                if (cah < gridx(1)) % x' smaller than min on grid
                    
                    disp('how can this be?')
                    
                end
                
                if (cah > gridx(nx)) % x' bigger than min on grid
                    
                    % set x' = max on grid
                    vp(ycc) = vpfun(jc+1,ycc,nx);
                    disp('you had more CAH than possible in the predefined grid--> adapt upper limit!')
                
                else
                
                    % interpolate V_x' at x' 
                    vp(ycc) = interp1(squeeze(gridx),squeeze(vpfun(jc+1,ycc,:)), cah);
                
                end
            end
            
            % Euler equation: RHS
            expvp=betta*sr(jc)*(1.0+r)*sum(pi(yc,:)*vp(:));
            
            % consumption
            cfun(jc,yc,xc)=invut(expvp);
            
            % check if budget constraint binds 
            if (cfun(jc,yc,xc)>=gridx(xc))
                
                cfun(jc,yc,xc) = gridx(xc);
                disp('shuu, borrowing constraint binds again')
                
            end
            
        end
        
        % income in current period/age (given shock yc)
        inc=epsi(jc)*netw*gridy(yc)+(1-epsi(jc))*pens;
        
        % assets at all xc:
        gridass(jc,yc,:)=(gridx-inc)/(1+r);
        
        % Update vfun and vpfun
        vpfun(jc,yc,:) = MUc(cfun(jc,yc,:));
        
        % Calculate implied value function tomorrow 
        for xc=1:nx
            
            v=zeros(ny,1);
            
            for ycc=1:ny
                
                % income tomorrow:
                incp1=epsi(jc+1)*netw*gridy(ycc)+(1-epsi(jc+1))*pens;
                
                % cah tomorrow
                cah=max(sqrt(eps),incp1+(1.0+r)*gridass(jc+1,ycc,xc));
                
                % this should never be the case:
                
                if ((cah + sqrt(eps)) < gridx(1))
                    
                    warning('How can this be ?');
                    
                end
                
                % linear interpolation:
                v(ycc)=func_intp(squeeze(gridx), squeeze(vfun(jc+1,ycc,:)), cah);
            
            end    
            
            % update value function with optimal consumption
            expv=sum(pi(yc,:)*v(:));
            vfun(jc,yc,xc)=U(cfun(jc,yc,xc))+betta*sr(jc)*expv;
            
        end    % end for xc
        
    end    % end for yc
    
end    % end for jc


% ---------------------------------------------------------------------
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
% ---------------------------------------------------------------------


% ---------------------------------------------------------------------
    function y=func_extrapol(x1,x2,y1,y2,x)
        
        % simple linear extrapolation
        
        m = (y2-y1)/(x2-x1);
        y = y1 + m*(x-x1);
        
    end
% ---------------------------------------------------------------------

end     % end function func_hh
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [Phi,PhiAss,ass]=func_aggr(gridx,gridsav,cfun,gridass)

global r nj nx ny pi gridy netw pens sr epsi pini frac totpop

disp('aggregation and cross-sectional measure');

% Compute Cross sectional distributions and aggregate variables
Phi = zeros(nj,ny,nx);          % distribution of assets conditional by age and shock
PhiAss = zeros(nx,1);             % distribution of assets

% Distribution of newborns over cash at hand
for yc=1:ny
    
    % income (wages and pensions) in current period/age:
    inc=epsi(1)*netw*gridy(yc)+(1-epsi(1))*pens;
    
    % initial cash-on-hand:
    cahini=inc;
    
    [vals,inds]=basefun(gridx,cahini,nx);
    Phi(1,yc,inds(1))=vals(1)*pini(yc)*frac(1);
    Phi(1,yc,inds(2))=vals(2)*pini(yc)*frac(1);
end;

for jc=2:nj
    TT = zeros(ny,nx,ny,nx);    % transfer function
    
    for xc=1:nx
        for yc=1:ny
            for ycc=1:ny
                
                % income (wages and pensions) in current period/age:
                inc=epsi(jc)*netw*gridy(ycc)+(1-epsi(jc))*pens;
                
                % cash on hand: x=a*(1+r)+y = s(-1)*(1+r)+y;
                %cah=inc+(1.0+r)*gridsav(xc);
                cah=inc+(1.0+r)*gridass(xc);
                
                [vals,inds]=basefun(gridx, cah,nx);
                
                TT(ycc,inds(1),yc,xc)=vals(1)*pi(yc,ycc);
                TT(ycc,inds(2),yc,xc)=vals(2)*pi(yc,ycc);
            end;    
        end;    
    end;    
    
    for xc=1:nx
        for yc=1:ny
            for xcc=1:nx
                for ycc=1:ny
                    % transfer distribution:
                    Phi(jc,ycc,xcc)=Phi(jc,ycc,xcc)+Phi(jc-1,yc,xc)*TT(ycc,xcc,yc,xc)*sr(jc-1);
                end;
            end;
        end;
    end;
    
end;    % end for jc


% Check if distribution sums to 1
sumprob=sum(sum(sum(Phi(:,:,:))));
if ( ( sumprob < 0.999 ) || ( sumprob > 1.001) ),
    beep; beep; beep;
    warning('distribution fucked up');
end;

% Check if Grid is Big enough
sumprob=sum(sum(Phi(:,:,nx)));
if (sumprob > 0.001 ),
    beep; beep; beep;
    warning('grid too small -- increase your grid');
    pause
end;

ass=0.0;
cons=0.0;
lab=0.0;
ret=0.0;

% aggregation
for jc=1:nj
    for yc=1:ny
        for xc=1:nx,
            PhiAss(xc)=PhiAss(xc)+Phi(jc,yc,xc);
            
            % asset holdings = capital stock in general equilibrium
            ass=ass+totpop*Phi(jc,yc,xc)*gridass(xc);
            
            cons=cons+totpop*Phi(jc,yc,xc)*cfun(jc,yc,xc);
            
            lab=lab+totpop*Phi(jc,yc,xc)*gridy(yc)*epsi(jc);
            ret=ret+totpop*Phi(jc,yc,xc)*gridy(yc)*(1.0-epsi(jc));
        end;
    end;
end;


% ---------------------------------------------------------------------
    function [vals,inds]=basefun(grid_x,x,nx)
        % this subroutine returns the values and the indices of the two basis
        % functions that are positive on a given x in the grid_x
        
        % MF function to lookup the current position
        i=lookup(grid_x,x,0);
        
        if ( (i+1)> nx-1),
            inds(1)=nx;
            inds(2)=nx;
            vals(2)=0.0;
            vals(1)=1.0;
        elseif (i==0),
            inds(1)=1;
            inds(2)=1;
            vals(1)=1.0;
            vals(2)=0.0;
        else
            inds(1)=i;
            inds(2)=i+1;
            dist = grid_x(i+1)-grid_x(i);
            vals(2)=( x-grid_x(i) )/dist;
            vals(1)=( grid_x(i+1)-x )/dist;
        end;
        
    end 	% end function basefun
% ---------------------------------------------------------------------




end     % end function func_aggr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [labinclife,inclife,asslife,conslife,vallife] = lcprofile(Phi,gridass,cfun,vfun)

global r nj nx ny gridy netw pens epsi frac

disp('life-cycle profiles')

asslife=zeros(nj,1);
inclife=zeros(nj,1);
labinclife=zeros(nj,1);
conslife=zeros(nj,1);
vallife=zeros(nj,1);

for jc=1:nj,
    for yc=1:ny
        for xc=1:nx,
            asslife(jc)=asslife(jc)+Phi(jc,yc,xc)*gridass(jc,yc,xc)/frac(jc);
            conslife(jc)=conslife(jc)+Phi(jc,yc,xc)*cfun(jc,yc,xc)/frac(jc);
            
            inc=epsi(jc)*netw.*gridy(yc)+(1-epsi(jc))*pens;
            labinclife(jc)=labinclife(jc)+Phi(jc,yc,xc)*inc/frac(jc);
            inclife(jc)=inclife(jc)+Phi(jc,yc,xc)*(r*gridass(jc,yc,xc)+inc)/frac(jc);
            
            vallife(jc)=vallife(jc)+Phi(jc,yc,xc)*vfun(jc,yc,xc)/frac(jc);
        end;
    end;
end;

end     % end function lcprofile
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

end  % end function mchain
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function u = U(c)
global tetta

% utility
if (abs(tetta-1-0)<sqrt(eps)),
    u = log(c);
else
    u = c.^(1.0-tetta)/(1.0-tetta);
end;
end     % end function U
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function muc=MUc(c)
global tetta

% maringal utility
muc = c.^(-tetta);
end     % end function MUc
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function invut=invut(marg)
global tetta

% invert utility for c
invut=marg.^(-1.0/tetta);
end     % end function invut
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function grd = makegrid(x1,x2,n,c)
% makes curved grid according to curvature parameter c
scale=x2-x1;
grd(1)=x1;
grd(n)=x2;
for i=2:n-1
    grd(i)=x1+scale*((i-1.0)/(n-1.0))^c;
end;
end     % end function makegrid
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
