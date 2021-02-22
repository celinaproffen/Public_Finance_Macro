function [g]=ginindex(a)
%  PURPUSE: Compute the Gini index for a given vector of income a and plot 
%           the correspnging Lorentz curve.                           
%--------------------------------------------------------------------------
%  USAGE:   [g]=ginindex(a)
%  Where:                                                                 
%            g = The Gini index.      
%            a = Vector of the observation to study      
%--------------------------------------------------------------------------
% REFERENCES: Portnov and Felsenstein (not dated):"Measures of Regional 
%             Inequality for Small Countries"
%--------------------------------------------------------------------------
% AUTHOR:  
% Karim Aroussi, 
% november 18, 2006.                                                  
% aroussi@em-lyon.com

n=length(a);
mu=mean(a);
for i=1:n
    if (a(i,1)<0)    error('Income value can not be negatif. Please verify your data'); end; 
end
%------------ calculate the Gini index ------------------*
s=0;
for i=1:n
    for j=1:n
        s=s+abs(a(i,1)-a(j,1));
    end
end
d=2*mu*((n^2));
g=s/d;
%-------------- plot the Lorentz curve -------------------*
k=sort(a);
kmax=max(k);
somme=sum(a);
z=zeros(n+1,1);
p=zeros(n+1,1);
for i=1:n
    z(i+1,1)= z(i,1)+(k(i,1)/(sum(a)));
    p(i+1,1)=p(i,1)+(1/(n));
end
x=[0 ;0.5; 1]; y=[0 ;0.5; 1]; % data created to plot the first bissectrice 
plot(p,z)
axis([0 1 0 1])
hold on;
plot(x,y,'red')
xlabel('Prcentage of population','FontSize',12)
ylabel('Percentage share of income','FontSize',12)
grid on
h = legend('Lorentz curve','Line of perfect equality',2); 

