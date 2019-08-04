function runNoisy

% this script runs the simulation of the AB model in landscapes
% with additive noise

addpath([pwd '/functions']);

%% parameters
gamma = [0.0500    0.0607    0.0737    0.0895    0.1087    0.1320    0.1603    0.1946];
%logspace(log10(0.05),log10(2),20); %  memory
mu  = [0.0005 0.001 0.005 0.01 0.025 0.05 0.5   1 ];
num = [    50    50    50  100   100  100 100 200 ];

beta = 5;             % amplification (keep beta*v0 < 0.2)
xi = 0.01;            % noise strength
T = 4000;             % simulation time  
N = 1000;             % number of cells
L = 20;               % domain
v = 0.02;             % swim velocity
saveInt = T/4;        % every how many seconds to take a snapshot
dx = 5e-5;            % mesh size (dt*v0, OK for mu > 0.0005)
dt = 5e-4;            % time step (OK for gamma > 0.1) 

folder = '/disk2/Adam/chemotaxis/beta5';
profiles = @(mu,xi) ['prof_OU_corr_',num2str(mu),'_xi_',num2str(xi),'.mat'];

% initial condition (*)
x0 = randn(N,1)*L/100 + L/3;

%% run AB model

addpath('/disk2/Adam/chemotaxis/profiles') 
for l = 1:length(mu)     
    data = load(profiles(mu(l),xi));
    S = data.S(:,1:num(l));
    tmpNum = num(l);
    tmpFolder = [folder '/mu_' num2str(mu(l))];
    
    parpool(40)
    parfor k = 1:tmpNum%i = 1:num*length(gamma)
        tmpS = S(:,k)';
        for j = 1:length(gamma)
           
        % this is to have the parfors different seeds
        pause(0.2)
        time = clock;
        seed = ceil(time(6)*100);
        rng(seed) 
            
        par = {T,L,dt,dx,v,beta/gamma(j)^2,gamma(j),saveInt,N,tmpFolder,k,seed};
        funAB(tmpS,x0,par);
        end
    end
    delete(gcp)
end