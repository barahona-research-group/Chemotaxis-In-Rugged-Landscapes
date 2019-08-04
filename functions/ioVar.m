function ioVar

%% Parameters
dt     = 5e-4;   %discretisation step
gamma  = logspace(log10(0.001),log10(1000),10); %memory range of interest
mu     = 10;%[0.0005 10]; %spatial correlation (extreme values)
v0     = 0.02;
mu     = mu/v0;  %effective (temporal) correlation length
T      = 3000;   %length of trajectories
num    = 100;    %number of realisations

% Note: here the length of trajectory is num x T. In the AB simulations we
% also need to capture good stats of the run-tumble dynamics!

%irrelevant here
sigma  = 1;      %noise strength
beta   = 1;
slope  = 0;      %slope
seed   = 0;      %random seed

avg = []; se = [];
for l = 1:length(mu)
    l
    
    % OU profiles
    par = {T,slope,sigma,mu(l),dt,num,seed,0,0};   
    S = profileOU(par);
    S = S(:,2:end);
    
    m = zeros(1,3);
    Delta = zeros(1,T/dt*size(S,2));
    meanDelta = zeros(1,length(gamma));
    stdDelta = zeros(1,length(gamma));

    for k = 1:length(gamma)      
        for j = 1:size(S,2)
            for i = 1:T/dt
                [m, Lambda] = Kt(m,S(i,j),gamma(k),dt);                
                Delta(T/dt*(j-1) + i) = beta/gamma(k)*Lambda;
            end
        end

        x = Delta(end/2-1:end); %assume this is stationary
        x = abs(x).^2; 
        meanDelta(k) = mean(x);
        stdDelta(k) = std(x);
    end

    avg = [avg; meanDelta]; se = [se; stdDelta];
end

%% exact
G = 0:0.1:5;
sigLam = G.*(G+3)./(8*(G+1).^3);
plot(G,sigLam)
hold on

%% plot
for k = 1:length(mu)
    fac = sigma;
    plot(gamma/mu(k),avg(k,:)/fac,'o')
end
xlabel('\Gamma')
ylabel('\sigma_\Lambda(\Gamma)')
xlim([0 5])