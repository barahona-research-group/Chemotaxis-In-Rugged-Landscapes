function procDet

% this script runs the simulation of the AB and KS models in deterministic
% linear landscape

addpath([pwd '/functions']);

%% common parameters
gamma = logspace(log10(0.05),log10(2),20);  % memory (for test try [0.5 1 2])
T = 3000;             %foraging time
v = 0.02;             %swim velocity
N = 50000;            % number of cells
folder = '/disk2/Adam/chemotaxis/Det/beta_20';
deter = 1;            % 1 : deterministic, 0 : noisy
Nprof = 0;            % number of noisy profiles
L = 20;               % length of domain
saveInt = T/3;        % every how many seconds to take a snapshot

% compute drift velocity, max Lambda for each gamma
figure
par = {T,L,v,gamma,N,folder,deter,Nprof};
procAB_driftVel(par)

% compute sample paths for a specific gamma
figure
par = {T,L,v,gamma(10),N,folder,deter,saveInt};
procAB_samplePaths(par)