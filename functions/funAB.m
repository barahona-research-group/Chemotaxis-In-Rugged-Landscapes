function [X,snaps,T,N,beta,gamma] = funAB(lig,x0,par)
       
%% Parameters
T = par{1};   %foraging time
L = par{2};   %domain 
dt = par{3};   %timegrid
dx = par{4};   %space grid
v = par{5};  % swimming velocity
beta = par{6};
gamma = par{7};
saveInt = par{8};
N = par{9};
folder = par{10};
profNum = par{11};
seed = par{12};

%% initial conditions
dir = 2*(rand(N,1)>0.5)-1; %randomised initial direction
x = x0; %initial positions
S0 = lig(ceil(x0/dx));
m = zeros(N,3);
m(:,1) = S0*gamma;
m(:,2) = S0*gamma^2;
m(:,3) = 2*S0*gamma^3;

%% check if file already exists
filename = ['/AB','_T_',num2str(T),'_N_',num2str(N),'_gamma_',num2str(gamma),'_',num2str(profNum),'.mat'];
if exist(folder,'file') ~= 7; mkdir(folder); end
if exist([folder filename],'file') == 2; disp('OK'); return; end

%% Run
fprintf('\n Counter: 0 \n');
snaps = zeros(1,T/saveInt+1);                % time vector of snapshots
X = zeros(N,T/saveInt+1); X(:,1) = x;        % cell positions
Dir = zeros(N,T/saveInt+1); Dir(:,1) = dir;  % run directions    
Deltat = zeros(N,T/saveInt);                 % tumbling response
% St = zeros(N,T/saveInt+1); St(:,1) = S0;   % Attractant concentrations

for t = 1:T/dt
    S = lig(ceil(x/dx))'; %ligand concentration at cell positions
    
    %tumbling response
    [m, Delta] = Kt(m,S,gamma,dt);
    lambda = max(0,1 - beta*Delta);   
    p = dt * lambda; % tumbling prob in interval of length dt

    % with probability p, tumble, and then update position
    dir = dir.*sign(rand(N,1) - p);
    x = x + dir*v*dt;

    % boundary conditions
    x = abs(x); %reflect from 0
    dir = dir.*sign(x);
    x = L - abs(L - x); %reflect from L
    dir = dir.*sign(L - x);  

    % collect data
    if  mod(t, saveInt/dt) == 0
        fprintf('\n Counter: %4.0f \n',t*dt)       
        snaps(t/saveInt*dt+1) = t*dt; 
        X(:,t/saveInt*dt+1) = x;  
        Deltat(:,t/saveInt*dt) = Delta;
        Dir(:,t/saveInt*dt+1) = dir;
%         St() = [St,S];
    end  
end 

%% save data
save([folder filename],'X','Deltat','snaps','T','N','beta','gamma','seed','-v7.3')