function rho = funKS(S,rho0,par)

%% Parameters
T = par{1};   %foraging time
L = par{2};   %domain 
dt = par{3};   %timegrid
dx = par{4};   %space grid
v = par{5};
beta = par{6};
gamma = par{7};
saveInt = par{8};
folder = par{9};
profNum = par{10};

% coarse-grained parameters
nu = v^2/2;  %diffusion
chi = 2*beta*gamma^2*nu/( 1 + 2*gamma )^3; %advection

M = L/dx - 1; N = T/dt + 1;
x = linspace(0,L,M+2);
S = interp1(linspace(0,L,length(S)),S,x);
q = dt/dx; p = nu*dt/dx^2;

% Diffusion condition
cfl = nu*dt/dx^2;
disp(['Max dt is ',num2str(dx^2/(2*nu)), ' with current M.'])
if cfl >= 1/2
    error('CFL condition is not satisfied! Decrease dt or increase h')
end

%% Matrices
% Difference operator
D = sparse(1:M-1,1:M-1,ones(1,M-1),M,M-1) ...
  - sparse(2:M,  1:M-1,ones(1,M-1),M,M-1);

% Tri-diagonal matrix for Laplace operator with Neumann BC
Lap = p*( diag(ones(M-1,1),-1) - 2*diag(ones(M,1),0) + diag(ones(M-1,1),1) );
A = eye(M) + Lap;
A(1,1) = 1-p; A(end,end) = 1-p;

%% Compute KS evolution
dS = diff(S)/dx;
v = chi*repmat(dS,1,N-1);
rho = zeros(M,N); rho(:,1) = rho0;
snaps = 0;
Rho = rho0;     
for t = 1:N-1        
    bx = sparse(1:M-1,2:M,v(:,t) < 0,M-1,M) ...
        + sparse(1:M-1,1:M-1,v(:,t) >= 0,M-1,M);
    B = q*D*diag(v(:,t))*bx;
        
    % compute (m^k+1)^i+1 using (v^k+1)^i 
    rho(:,t+1) = (A - B)*rho(:,t);
    
    % collect data
    if  mod(t, saveInt/dt) == 0       
        snaps = [snaps, t/dt]; 
        Rho = [Rho, rho(:,t+1)];
    end 
end

%% save data
if exist(folder,'file') ~= 7; mkdir(folder); end
filename = ['/KS','_T_',num2str(T),'_gamma_',num2str(gamma),'_',num2str(profNum),'.mat'];
save([folder filename],'Rho','snaps','T','beta','gamma','-v7.3')