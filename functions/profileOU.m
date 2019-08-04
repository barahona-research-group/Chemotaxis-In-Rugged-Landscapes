function S = profileOU(par)

% Generate a number of realisations of the linear gradient with additive,
% spatially correlated noise. 
%
% Noise is simulated based on the Ornstein-Uhlenbeck (OU) process
%
%   dXt = (aXt + b) + xi dWt, 
%
% where W is a unit Wiener process. Here  b = 0; a = -1/mu; and
% xi = sqrt(2*sigma/mu) in order to normalise the noise to have the 
% same amplitude.

%% parameters
L      = par{1}; %domain length
slope  = par{2}; %slope
sigma  = par{3}; %noise strength
mu     = par{4}; %correlation length [mm], 0.01-0.1mm=100mum
dx     = par{5}; %discretisation step
num    = par{6}; %number of realisations
seed   = par{7}; %random seed
folder = par{8};
output = par{9};
a = -1/mu; 
b  = 0;          
xi = sqrt(2*sigma/mu); %same amplitude
% xi = sqrt(2*sigma/mu^2); %same psd

%% deterministic part
Sfun = @(x) slope * x;
S = repmat(Sfun((-L/2:dx:L/2-dx)'),1,num+1);

%% Simulate OU process 

% using SDE toolbox
% root = '/Users/adamgosztolai/Documents/Research/Software/SDE_Toolbox_1.4.1';
% root = '/home/adam/toolboxes/SDE_Toolbox_1.4.1';
% addpath(root);
% addpath([root,'/models_library']);
% x0 = 0;       % the SDE initial condition
% numVars = 1;  % the dimension of the SDE
% xhat = SDE_euler([x0,a,b,xi],'M2',x0:dx:L-dx,numVars,num,'Ito',seed);

%alternatively 
xhat = zeros(L/dx,num);
for j = 1:num
    rng(j)  
    for i = 1:L/dx-1
        xhat(i+1,j) = xhat(i,j) + a*xhat(i,j)*dx + xi*sqrt(dx)*randn;
    end
end

%% deterministic + additive OU noise
S(:,2:end) = S(:,2:end) + xhat;

%% save
if output == 1
    if exist(folder,'file') ~= 7; mkdir(folder); end
    filename = ['/prof_OU_corr_',num2str(mu),'_xi_',num2str(sigma),'.mat'];
    save([folder filename],'S','L','sigma','mu','dx','num','slope')
end