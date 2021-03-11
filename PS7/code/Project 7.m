global alpha beta tau lambda T N s Phi phi eta varrho
alpha = 0.3;
beta = 0.99^40;
tau = 0;
lambda=0.5;
T =50000;
N=1000; % number of simulations
s=zeros(1,T)
Phi=zeros(1,T)
phi=zeros(N,T)
eta=zeros(N,T)
varrho=zeros(N,T)

function phi=phi(alpha,lambda,varrho,eta,tau)

function s= saving(beta,phi)
    s=(beta*phi)/(1+beta*phi)
end
 % saving function 
function s = phi(beta, alpha, tau, eta, lambda)
        for n=1:N
            phi(n,t)=1/(1+((1-alpha)*(lambda*eta(n,t+1)/(alpha*(1+lambda)*varrho(n,t+1))))
        end
        Phi(1,t)=sum(phi(t,n))/N
end
end


        
            
        



