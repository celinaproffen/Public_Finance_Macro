% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% PROJECT 6: Exercise 1
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function exercise

close all

global alpha delta

delta = 0.05; % orig = 0.05
alpha = 0.33; % 0.33; % 0.1;
ret = 0.02; % initial guess

tol = 1e-4;
maxit = 10;
df = 0.1;

retold = ret;

% Remark: We deleted the graphs as we don't need them for the project.

% solution of OLG model:
ret = retold;

for it=1:maxit,
    [fval, ass, Y, wage_scale] = func_olg(ret);
    if (abs(fval)<tol),
        disp('convergence');
        break;
    else,
        ret = ret - df * fval;
        disp(['iteration #', num2str(it)]);
        disp(['guess of rate of return: ', num2str(ret)]);
        disp(['wage scaling factor: ', num2str(wage_scale)]);
    end;
end;
if (it>=maxit),
    warning('no convergence');
end;
disp(['equilibrium interest rate: ', num2str(ret)]);
%disp(['equilibrium contribution rate: ', num2str(tau)]);
disp(['equilibrium capital output ratio: ', num2str(ass/Y)]);

end

function [fval,ass,Y, wage_scale] = func_olg(ret)
   
    %%% NEW %%%
    global alpha delta

    mpk = ret + delta; % marginal product of capital

    wage_scale = func_firm(mpk); % marginal product of labor that depends on capital

    % idea: feed the wage level and current guess for rate of return
    % into the towards_olg function 

    ass = towards_olg(wage_scale,ret)

    [mpk,Y] = func_mpk_altern(ass); 
    %%% end NEW %%%
    
    retnew = mpk - delta 
    
    fval = ret - retnew; % Change in rate of return 
end


function [aggregate_cap_supply] = towards_olg(wage_scaling_factor,ret)

%close all

global nj ny

tic

opt_det=false;          % 1=deterministic model (1=true, 0=false in matlab)
opt_nosr=1;             % 1=no survival risk (if we want survival risk later, set this to 0)
opt_ny = 2;             % 1=Markov chain with number of states, ny=5, % 2=Markov chain with ny=2 (Krueger-Ludwig calibration)

% -------------------------------------------------------------------------
% SOLUTION

% calibration
func_calibr(opt_det,opt_nosr,opt_ny,wage_scaling_factor,ret);

% solution of household model
[gridx,gridsav,gridass,cfun,vfun] = func_hh;

% aggregation
[Phi,PhiAss,ass] = func_aggr(gridx,gridsav,cfun,gridass);

% average life-cycle profiles
[labinclife,inclife,asslife,conslife,vallife] = lcprofile(Phi,gridass,cfun,vfun);

%%% NEW %%%
aggregate_cap_supply = ass;
%%% end NEW %%%

% -------------------------------------------------------------------------

% We deleted the plots here as well as the tic-toc, as it was not needed
% for the project.

end     % end function func_main
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function func_calibr(opt_det,opt_nosr,opt_ny,wage_scaling_factor,ret)

global betta tetta nj jr nx ny pi gridy netw pens sr epsi curv pini frac pop totpop grdfac r

close all

%%% NEW %%%
r = ret; % now this is determined endogenously
%%% end NEW %%%
rho = 0.04;
betta = 1/(1+rho);
tetta = 2;

nj=80; % number of year a person lives
jr=45; % after 45 year people retire

nx=30;         % # of grid-points
curv=3.0;       % curvature of grid
grdfac=40;      % scaling factor of saving grid


% deterministic income component:
netw=1.0;
pens=0.4; % originally = 0.4;

% net = 1 --> post_government(yhh5) vs. netincome == 2 --> pre-government(yhh6)
netincome = 2; % 1

if (netincome == 1),
    columnhelp = 9;
elseif (netincome==2),
    columnhelp = 6;
end;
    
epsi=ones(nj,1);

inc = readmatrix('predicted_values.xlsx'); % we input income from predicted values file as a matrix in matlab.
epsi = inc(:,columnhelp); % we replace the normalized vector of one with the vector of our estimates (for pre or post gov repectively) 

%%% NEW %%%
epsi = epsi*wage_scaling_factor;
%%% end NEW %%%

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
        
        % ii) exact by eigenvector decomposition:
        [v,d] = eig(pi');
        v = v./sum(v);
        pini = v(:,1);        
        % take exponent and rescale income shocks such that mean is one:
        gridy=exp(gridy)';
        gridy = gridy / sum(gridy.*pini);
        
    else        
        % Alternative -- taken from Kr??ger and Ludwig (2007) using only two
        % states
        ny = 2;
        
        vary=0.08;    % taken from Storesletten, Telmer, Yaron
        
        % transition probability and variance given pre- or post gov income
        if (netincome==1),
            rhoeta=0.9049 ; % We replace with our estimate of transition probability from problem 2 
            epsil=0.7187; % shock for post-government income
        elseif (netincome ==2),
            rhoeta=0.9485; % We replace with our estimate of transition probability from problem 2 
            epsil=0.8424; % shock for pre-government income
        end;
        
        % Markov chain
        [pini,pi,gridy]=mchain(rhoeta,epsil);
    end;
    
end;

end     % end function func_calibr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [gridx,gridsav,gridass,cfun,vfun] = func_hh

global betta tetta r nj nx ny pi gridy netw pens sr epsi curv grdfac

disp('solution of household model');

% grids and decisions rules:
gridx = zeros(nj,ny,nx); % years of life * possible distinct numbers for permanent shock realizations * possible cash-at-hand-states (three dimensional!)
gridsav = zeros(nx,1); % savings grid which we need in the alternative specification implemented here (one dimensional!)
gridass = zeros(nj,ny,nx); % to store associated assets => years of life * possible distinct numbers for permanent shock realizations * possible cash-at-hand-states (three dimensional!)
cfun = zeros(nj,ny,nx); % policy function in all states and ages for different cash-at-hand levels
vfun = zeros(nj,ny,nx); % value function in all states and ages for different cash-at-hand levels
vpfun = zeros(nx,1); % first derivative of the value function (w.r.t.consumption) (? why is the formula 8.62 in the dynmaic macro notes different?)
vptrans = zeros(nj,ny,nx); % transformed value function in all states and ages for different cash-at-hand levels

% savings grid: hold it constant:
maxsav=grdfac;
gridsav(2:nx)=makegrid(0.0,grdfac,nx-1,curv);
gridsav(1)=0.0;

% for all income states define possible cash-at-hand holdings in a 3
% dimensional grid (age, state, amount).
% Do the same for savings asset holdings, consumption functions, value
% functions, transformed value functions, etc.
for yc=1:ny
    % cash-on-hand grid at nj; thus is an nj*1 vector.
    inc = epsi(nj)*netw*gridy(yc)+(1-epsi(nj))*pens;
    % when you are young nj<jr your receive income 1 multiplies by the
    % shock; we don't care about transition probabilities yet but will
    % implement them later!
    
    % in case of no pension system (which is not the case in Exercise 1, but could happen), assume some minimum cash on hand:
    % This is needed to explain precautionary motives in t-1
    minx=max(inc,sqrt(eps));
    maxx=gridsav(nx)*(1.0+r)+inc;
    gridx(nj,yc,:)=linspace(minx,maxx,nx); % This is really important: 
    % You basically built up a 3 dimensional grid, for all ages and states of shocks you construct a gridspace for cash-at-hand
    % The minimum cash at hand you can have is somewhat above zero as we
    % are using the theory of INADA utility functions
    
    % Final period Consumption function, asset holdings (before interest rates realize, see instructions in project), value function, including derivative
    cfun(nj,yc,:)=gridx(nj,yc,:);
    gridass(nj,yc,:)=(gridx(nj,yc,:)-inc)/(1+r);
    vfun(nj,yc,:)=U(cfun(nj,yc,:));
    vpfun(:)=MUc(cfun(nj,yc,:));
    vptrans(nj,yc,:)=vpfun.^(-1.0/tetta);
end;

% Iterate Backwards
for jc=nj-1:-1:1,
    
    for yc=1:ny % state from which you are coming/ shocks today
        
        for xc=2:nx, % nx is the number of grid points of savings
            vp=zeros(2,1); % here it is two as we only have two income states.
            
            for ycc=1:ny,
                % income tomorrow:
                incp1=epsi(jc+1)*netw*gridy(ycc)+(1-epsi(jc+1))*pens;
                
                % Maximum cash on hand tomorrow:
                % in case of zero savings and no pension system assume some
                % minimum cash on hand
                cah=max(sqrt(eps),incp1+(1.0+r)*gridsav(xc));
                
                % Interpolate derivative of value function
                if ( cah<gridx(jc+1,ycc,1)),
                    disp('how can this be?') % if this applies you have calculated that the person 
                    %receives LESS than the minimum value on the grid of possible CAH for ages and
                    % states of income shocks.. This would make no sense
                end;
                if ( cah>gridx(jc+1,ycc,nx) ),
                    % if out of bounds simply set it to decision at nx:
                    vptr = vptrans(jc+1,ycc,nx);
                else
                    vptr = interp1(squeeze(gridx(jc+1,ycc,:)),squeeze(vptrans(jc+1,ycc,:)),cah);
                    % this gives you the array of values in the transformed
                    % value function, that would arise if your cash at hand
                    % cah was in statae ycc and savings were at any of the
                    % possible grip points gridsav(xc)
                end;
                vp(ycc)=vptr.^(-tetta); 
                % gives you the actual vaue that results from a given state
                % and savings decision
            end;
            
            % Euler equation: RHS
            % expected consumption next period (discounted and adjusted for
            % time prefeences
            expvp=betta*sr(jc)*(1.0+r)*sum(pi(yc,:)*vp(:));
            
            % consumption (which was a function of the marginal utility of
            % consumption tomorrow, which we need to revert back to actual
            % consumption values here)
            % !!!!!
            cfun(jc,yc,xc)=invut(expvp);
            
            % endogenous x-grid for cah:
            gridx(jc,yc,xc)=gridsav(xc)+cfun(jc,yc,xc);
        end;
        
        % income (wages and pensions) in current period/age:
        inc=epsi(jc)*netw*gridy(yc)+(1-epsi(jc))*pens;
        
        % decision at minx
        % notice: correction required for welfare calculation
        % the above is actually slightly inefficient because xmin
        % can be explicitly computed, then gridsav would be age and
        % state dependent.
        % Basically we adjust the lowest grid point for x (CAH) here
        % ASK! Why do we need this?
        minx=max(inc,sqrt(eps));
        if (minx<gridx(jc,yc,2)),
            gridx(jc,yc,1)=minx;
        else    % set it to some arbitrary fracion of x(2)
            gridx(jc,yc,1)=0.9*gridx(jc,yc,2);
        end;
        
        % Compute optimal consumption for minx
        % The minimum consumption is set to the never reached minimum level
        % of CAH that could ever be realized
        cfun(jc,yc,1)=gridx(jc,yc,1);
        
        % assets at all xc:
        gridass(jc,yc,:)=(gridx(jc,yc,:)-inc)/(1+r);
        
        % Update vfun and vpfun
        vpfun(:)=MUc(cfun(jc,yc,:));
        vptrans(jc,yc,:)=vpfun(:).^(-1.0/tetta);
        
        % Calculate value function
        for xc=1:nx,
            
            v=zeros(2,1);
            for ycc=1:ny,
                % income tomorrow:
                incp1=epsi(jc+1)*netw*gridy(ycc)+(1-epsi(jc+1))*pens;
                
                % cah tomorrow
                cah=max(sqrt(eps),incp1+(1.0+r)*gridsav(xc));
                
                % this should never be the case:
                if ((cah+0.0001)<gridx(jc+1,ycc,1)),
                    warning('How can this be ?');
                end;
                % linear interpolation:
                v(ycc)=func_intp(squeeze(gridx(jc+1,ycc,:)),squeeze(vfun(jc+1,ycc,:)),cah);
            end;    % end for ycc
            
            % update value function
            expv=sum(pi(yc,:)*v(:));
            % this is the bellmann equation; what we fill in here is the
            % optimal consumption for each level of CAH
            vfun(jc,yc,xc)=U(cfun(jc,yc,xc))+betta*sr(jc)*expv;
        end;    % end for xc
        
    end;    % end for yc
    
end;    % end for jc


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
    
    [vals,inds]=basefun(gridx(1,yc,:),cahini,nx);
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
                cah=inc+(1.0+r)*gridsav(xc);
                
                [vals,inds]=basefun(gridx(jc,ycc,:),cah,nx);
                
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
            ass=ass+totpop*Phi(jc,yc,xc)*gridsav(xc);
            
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
        % I think this gives you the closest level in grid_x to x (here
        % cah) in term of the index of grid_x where to find that cah-level
        i=lookup(grid_x,x,0);
        
        if ( (i+1)>nx),
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
        
    end 	% end function basefun: you got the values and indices of the grid_x which u needed!
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

% -------------------------------------------
function wage=func_firm(mpk)

global alpha delta

k = (alpha/mpk)^(1/(1-alpha));
wage = (1-alpha)*k^alpha;

end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ----------------------------------------
function [mpk,Y] = func_mpk_altern(ass)

global alpha jr pop

L_help = pop(1:jr);
L = sum(L_help);
Y = ass.^alpha * L.^(1-alpha);
ky = ass./Y;
mpk = alpha * ky.^(-1);

end