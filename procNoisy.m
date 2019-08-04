function procNoisy

% this script runs the simulation of the AB and KS models in landscapes
% with additive noise

addpath([pwd '/functions']);

%% common parameters
gamma = logspace(log10(0.05),log10(2),20); %  memory
%gamma = [0.0500    0.0607    0.0737    0.0895    0.1087    0.1320    0.1603    0.1946];
mu  = [0.0005 0.001  0.005  0.01  0.025   0.05   0.5     1 ];
num = [    50    50     50   100    100    100   100   200 ];
N   = [  1000  1000   1000  1000    1000  1000  1000  1000 ];

T = 4000;             % foraging time
v = 0.02;             % swim velocity
folder = '/disk2/Adam/chemotaxis/beta5_s0.01';
deter = 0;            % 1 : deterministic, 0 : noisy
L = 20;               % length of domain
saveInt = T/4;        % every how many seconds to take a snapshot

% compute drift velocity, max Lambda for each gamma
figure
for i = 1:length(mu)
    tmpFolder = [folder '/mu_' num2str(mu(i))];
    par = {T,L,v,gamma,N(i),tmpFolder,deter,num(i),mu(i)};
    procAB_driftVel(par);
end

% compute sample paths for a specific gamma
par = {T,L,v,gamma(10),N,folder,deter,saveInt};
procAB_samplePaths(par)

% compute distance to Gaussian
figure
for i = 1:length(mu)
    tmpFolder = [folder '/mu_' num2str(mu(i))];
    par = {T,L,v,gamma,N(i),tmpFolder,deter,num(i),mu(i)};
    procAB_distToGauss(par)
end