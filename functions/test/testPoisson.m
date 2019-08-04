function testPoisson
%This script implements an approximation of the homogeneous Poisson process
%X ~ Pois(lambda), valid for lambda << 1.
%
%we test the accuracy of the discretisation for different dt.

%parameters
dt = [0.0001 0.0005 0.001 0.005 ]; %discretisation
lambda = .1; %rate
T = 50000; %simulation time

%simulate Poisson process
for j = 1:length(dt)
    
    X = [];
    for i = 1:T/dt(j)
        p = dt(j) * lambda; %probability of event in [t, t+dt]
        X(i) = (rand - p) < 0;
    end

    waitingtimes = diff(dt(j)*find(X));

    %mean
    mu(j) = mean(waitingtimes); 
    
    %confidence intervals
    SE = std(waitingtimes)/sqrt(length(waitingtimes));  % Standard Error
    ts = tinv([0.025  0.975],length(waitingtimes)-1);   % T-Score
    CI(j,:) = SE;
end

%plot discretisation step against mean waiting time
shadedErrorBar(log(dt),mu,CI)
xlabel('log(dt)')
ylabel('mean waiting time (1/rate)')