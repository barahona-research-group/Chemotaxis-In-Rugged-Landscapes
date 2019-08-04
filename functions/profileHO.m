function S = profileHO(par)

% Generate a number of realisations of the linear gradient with additive,
% spatially correlated noise. 

% Noise is simulated based on the stochastic harmonic oscillator model.
%
%   M ddx1 + B dx1 + K x1 = xi Phi, 
%
% where Phi is the unit white noise process. Equivalently,
%
%   dx1 = x2; 
%   dx2 = - K/M x1 - B/M x2 + xi/M Phi
%
% where the structural parameters read 
%
%   B = 1; K = 1/mu; xi = sqrt(2*sigma/mu)
%
% in order to normalise the noise to have the same amplitude.

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
M      = par{10}; %smoothing (for HO)
B      = 1;
K      = 1/mu;
xi = sqrt(2*sigma/mu);

%% deterministic part
Sfun = @(x) slope * x;
S = repmat(Sfun((-L/2:dx:L/2-dx)'),1,num+1);

%% Simulate HO process

%using SDE toolbox
% root = '/Users/adamgosztolai/Documents/Research/Software/SDE_Toolbox_1.4.1';
% root = '/home/adam/toolboxes/SDE_Toolbox_1.4.1';
% addpath(root);
% addpath([root,'/models_library'])
% numVars = 2; % the dimension of the SDE
% x0 = [0 0]; % the SDE initial condition
% pars = [0,0,0,-1,K/M,B/M]; % (?1, ?2, ?11, ?12, ?21, ?22)
% xi = [0, xi/M]; % the SDE structural parameter ?sigma?
% xhat = SDE_euler([x0,pars,xi],'M9',0:dx:L-dx,numVars,num,'Ito',seed);
% xhat = xhat(:,1:2:end-1);
    
%alternatively 
xhat = zeros(L/dx,num);
yhat = zeros(L/dx,1);
for j = 1:num
    rng(j)   
    for i = 1:L/dx-1
        xhat(i+1,j) = xhat(i,j) + yhat(i)*dx;
        yhat(i+1)   = yhat(i)*(1 - B/M*dx) - K/M*xhat(i,j)*dx + xi*sqrt(dx)*randn;
    end
end

S(:,2:end) = S(:,2:end) + xhat;

%% save
if output == 1
    filename = ['/prof_HO_corr_',num2str(mu),...
        '_xi_',num2str(sigma),'_in_',num2str(M),'.mat'];
    save([folder filename],'S','L','sigma','mu','dx','num','slope','M')
end