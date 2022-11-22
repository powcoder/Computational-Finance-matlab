function [X,PNL]=deltahedge(S0,r,mu,sigma,T,n,N)
% S0 is the initial stock price, 
% K is strike price,  
% r is free interest rate, 
% T is time to maturity, 
% n is the number of time discretisation steps,
% N is the sample size
% mu is the asset's rate of return 
% sigma is the asset's volatility
%
strike = linspace(80,120,21);
nstep = power(2,n);
delt = T/nstep;
const = (mu - 0.5*sigma*sigma);
tstep = linspace(0,T,nstep+1);
tspace = ones(N,1)*tstep;
BMval = normrnd(0,1,[N,nstep+1]);
BMval(1:N,1) = 0;
Wt = cumsum(BMval,2);
ST = S0*exp(const*tspace+sigma*sqrt(delt)*Wt);
%Compute delta at each readjustment time
tt = T-tstep(1:nstep);
X = zeros(N,21);
PNL = zeros(N,21);
for i=1:N
    for j=1:21
        price = blsprice(S0,strike(j),r,T,sigma);
        d1 = (log(ST(i,1:nstep)./(strike(j)*exp(-r*tt))) + 0.5*sigma*sigma*tt)./(sigma*sqrt(tt));
        delta = normcdf(d1);
        temp = ST(i,2:nstep+1).*exp(-r*tstep(2:nstep+1)) - ST(i,1:nstep).*exp(-r*tstep(1:nstep));
        chng = delta.*temp;
        X(i,j) = exp(r*T)*(price +sum(chng));
        PNL(i,j) = X(i,j) - max(ST(i,nstep+1)-strike(j),0);
    end
end
    
    