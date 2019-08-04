function procAB_distToGauss(par)

%parameters
T = par{1};
gamma = par{4};
N = par{5};
folder = par{6};
Nprof = par{8};

range = 1:Nprof;
for k = 1:length(gamma) % loop over gamma
    for i = range % loop over number of samples
        disp(gamma(k))

        file = ['/AB_T_',num2str(T),'_N_',num2str(N),...
            '_gamma_',num2str(gamma(k)),'_',num2str(i),'.mat'];
        
        if exist([folder file], 'file') == 0 
            error(['No result file for ' folder file])
        else            
            %load data
            data = load([folder file]);
            X = data.X;
         
            % compute empirical drift speed (ignoring initial transient)
            [L2dist(i),~,sigma(i)] = distToGauss(X);                     
        end
    end
      
    distAvg(k) = mean(L2dist);
    sigmaAvg(k) = mean(sigma);
    
end

%% plot
subplot(2,1,1); plot(gamma,distAvg); hold on
subplot(2,1,2); plot(gamma,sigmaAvg); hold on

end

function [L2dist,mu,sigma] = distToGauss(X)
    %compute L2 distance between histogram X and continuous pdf Y
    [binValues, binEdges] = histcounts(X(1:5000,end),200,'Normalization','pdf');
	binCentres = ( binEdges(1:end-1) + binEdges(2:end) )/2;
    hfit = fitdist(X(1:5000,end),'normal');
% 	hfit = histfit(X(:,end));
    mu = hfit.mu;
    sigma = hfit.sigma;
	normDistr = normpdf(binCentres,mu,sigma);
	L2dist = norm(normDistr - binValues);
end