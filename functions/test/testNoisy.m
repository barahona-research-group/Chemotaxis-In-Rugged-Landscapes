function testNoisy

% this script runs the simulation of the AB and KS models in landscapes
% with additive noise

addpath('../');

%% common parameters
gamma = [0.05 2]; %  memory
mu = [0.0005 10];
beta = 5;             % amplification (keep beta*v0 < 0.2)
xi = 0.001;           % noise strength
T = 3000;             % simulation time  
N = [100 500 1000 5000 10000 50000]; % number of cells
L = 20;               % domain
v = 0.02;             % swim velocity
saveInt = T/3;        % every how many seconds to take a snapshot
dx = 1e-4;            % mesh size (dt*v0, OK for mu > 0.0005)
dt = 1e-3;            % time step (OK for gamma > 0.1)

% Note: T and N need to be large enough to capture the run-tumble dynamics 
% in a given realisation of the landscape. We fix T = 3000. Thus for for 
% longer correlation lengths we need longer paths or more cells to capture 
% the given landscape. Then, set number of landscape realisations high 
% enough to capture the statistics of the landscapes. 

folder = '/disk2/Adam/chemotaxis/test/';
profiles = @(mu,xi) ['prof_OU_corr_',num2str(mu),'_xi_',num2str(xi),'.mat'];


%% run AB model

addpath('/disk2/Adam/chemotaxis/profiles/test') 

for l = 1:length(mu)     
    data = load(profiles(mu(l),xi));
    S = data.S(:,2)';
    
    parpool(length(N)*length(gamma))
    parfor i = 1:length(N)*length(gamma)
        k = ceil(i/length(gamma)); %N
        j = i-(k-1)*length(gamma); %gamma
           
        rng(i) % this is to have the parfors different seeds
            
        % initial condition (*)
        x0 = randn(N(k),1)*L/200 + L/3;
        par = {T,L,dt,dx,v,beta/gamma(j)^2,gamma(j),saveInt,N(k),folder,1,i};
        funAB(S,x0,par);
    end
    delete(gcp)
end

rmpath('/disk2/Adam/chemotaxis/profiles/test')