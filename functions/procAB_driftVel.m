function [vmean, vstd] = procAB_driftVel(par)

%parameters
T = par{1};
v0 = par{3};
gamma = par{4};
N = par{5};
folder = par{6};
deter = par{7};
Nprof = par{8};
nbins = 200;
if deter == 0
    mu = par{9};
end

%% empirical drift velocity
if deter == 1; range = 1; elseif deter == 0; range = 1:Nprof; end
counts = zeros(length(gamma),nbins); edges = 0; num = 0;
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
            B(k) = data.beta;
            
            % compute empirical drift speed (ignoring initial transient)
            DeltaX = X(:,end)-X(:,1);
            vmean_temp(i) = mean(DeltaX)/T;
            vstd_temp(i)  =  std(DeltaX)/T;                      
        end
    end
      
    avgHist = sum(counts,2)/num;
    vmean(k) = mean(vmean_temp);
    vstd(k) = mean(vstd_temp);
    
    clearvars vmean_temp vstd_temp
end

%% plot
shadedErrorBar(gamma,vmeannorm./fac,vstd./fac,'ro')