function procProfiles

%% parameters
xi = 1e-3;
M = 1e-2;
mu = 1; %correlation length [mm], 0.01-0.1mm=100mum
num    = 1;     %number of realisations
folder = '/disk2/Adam/chemotaxis/profiles/test';
taumax = 10; 

%% Generate profiles 
for i = 1:length(mu)
    for j = 1
    
        % load profile
        if j == 1
            data = load([folder, '/prof_OU_corr_',num2str(mu(i)), ...
                '_xi_',num2str(xi),'.mat']);
        elseif j == 2
            data = load([folder, '/prof_HO_corr_',num2str(mu(i)), ...
                '_xi_',num2str(xi),'_in_',num2str(M),'.mat']);
        end
        
        S = data.S;
        dx = data.dx;
        xi = data.sigma;
        
        par = {xi,dx,num,taumax};
   
        % obtain autocorrelation function 
        R = autoCorr(S,par);
        x = 0:dx:taumax-dx; 
        
        % plot
        f = fit(x',R,'exp1');
        plot(f,x,R)          
        hold on    

        % parameters
        f
    end
end

%% plot
% figure
% plot(linspace(0,L,size(S,1)),S(:,2))