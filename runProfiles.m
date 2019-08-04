function runProfiles

addpath([pwd '/functions'])

%% parameters
L      = 20;     %domain length
slope  = 1;      %slope
sigma  = 1e-2;   %noise strength
mu = [0.00025 0.0005 0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 1 10]; %correlation length [mm], 0.01-0.1mm=100mum
dx     = 1e-4;   %discretisation step
num    = 200;    %number of realisations
seed   = 0;      %random seed
folder = '/disk2/Adam/chemotaxis/profiles';
output = 1;      %save?
M      = 1e-4;   %smoothing (for HO)

%% Generate profiles
for i = 1:length(mu)
    i
    par = {L,slope,sigma,mu(i),dx,num,seed,folder,output,M};
    
    % OU profiles
    profileOU(par);
    
    % HO profiles
    profileHO(par);
end