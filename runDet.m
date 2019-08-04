function runDet

% this script runs the simulation of the AB and KS models in deterministic
% linear landscape

addpath([pwd '/functions']);

%% common parameters
gamma = logspace(log10(0.05),log10(2),20);  % memory 
beta = 20;            %amplification
alpha = 1;            %slope of gradient
T = 3000;             %foraging time
L = 20;               %domain
v = 0.02;             %swim velocity
saveInt = T/3;        %every how many seconds to take a snapshot
folder = ['/disk2/Adam/chemotaxis/Det/beta_',num2str(beta)];

%% run ABM model

N = 50000;  % number of cells
dx = 0.001; % < 0.001 else too much round-off
dt = 0.005; % time step [sec]

% profile
S = alpha*(0:dx:L);

% initial condition (*)
x0 = randn(N,1)*L/200 + L/3;

% run AB model
parpool(length(gamma))
parfor i = 1:length(gamma)
% for i = 1:length(gamma)
    
    % this is to have the parfors different seeds
    pause(0.2)
    time = clock;
    seed = ceil(time(6)*100);
    rng(seed)
    
    par = {T,L,dt,dx,v,beta/gamma(i)^2,gamma(i),saveInt,N,folder,1,seed};
    funAB(S,x0,par);

end
delete(gcp)

%% run KS model
    
dt = 2;     %timegrid
dx = 0.05;  %space grid

% profile
S = alpha*(0:dx:L);

% initial condition (compare to *)
rho0 = normpdf(dx:dx:L-dx,L/3,L/50)';

% run KS
parpool(length(gamma))
parfor i = 1:length(gamma)
% for i = 1:length(gamma)
    
    par = {T,L,dt,dx,v,beta/gamma(i)^2,gamma(i),saveInt,folder,1};
    funKS(S,rho0,par);

end
delete(gcp)