len = 20; % number of samples to take
tau = [0.1 0.5 1 5];
T = 3000;

%AB solution (change)
fileAB = ['pde_T_',num2str(T),'_epsi_',num2str(tau(k)),'_',num2str(i),'.mat'];

% folders
folder ='/disk4/Adam/chemotaxis/Paper1/Fig4_data/';

for k = 1:length(tau)
    disp(tau(k))
    j = 0;
for i = 1:len  
    file = ['pde_T_',num2str(T),'_epsi_',num2str(tau(k)),'_',num2str(i),'.mat'];
    if exist([folder file], 'file') ~= 0 
        j = j+1;
        load([folder file])
        rhoTemp(:,j) = rho(:,end);
    else
        error('file does not exist')
    end

    clearvars -except i j k epsi mu len folder T tau file rhoTemp
end
    rhoAvg(:,k) = mean(rhoTemp,2);
    clear rhoTemp
    plot(linspace(0,10,length(rhoAvg(:,k))),rhoAvg(:,k))
    hold on
end

load(fileAB)
for k = 1:length(tau)
    
    L2dist(k) = distToAB(X,rhoAvg(:,k));
end

plot(tau,L2dist)

function L2dist = distToAB(X,rho,x)
    %compute L2 distance between histogram X and continuous pdf Y
    h = histogram(X(:,end),200,'Normalization','pdf');
	binCentres = ( h.BinEdges(1:end-1) + h.BinEdges(2:end) )/2;
	binValues = h.Values;
    
    rho = interp1(rho,x,binCentres);
	L2dist = norm(rho - binValues);
end