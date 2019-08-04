function procTestNoisy

% this script runs the simulation of the AB and KS models in deterministic
% linear landscape

addpath('../');

%% common parameters
gamma = [0.05 2]; % memory 
T = 3000;             %foraging time
v = 0.02;             %swim velocity
N = [100 500 1000 5000 10000 50000 100000]; % number of cells
folder = '/disk2/Adam/chemotaxis/test';
deter = 1;            % 1 : deterministic, 0 : noisy
Nprof = 1;            % number of noisy profiles
L = 20;

for i = 1:length(N)
    % compute drift velocity, max Lambda for each gamma
    par = {T,L,v,gamma,N(i),folder,deter,Nprof};
    [a, b] = procAB_driftVel(par);
    vmean(i,:) = a;
    vstd(i,:) = b;
end
plot(vmean)
