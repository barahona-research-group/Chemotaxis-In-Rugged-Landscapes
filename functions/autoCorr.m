function R = autoCorr(y,par)

%% parameters
sigma  = par{1};
dx     = par{2};
num    = par{3};
taumax = par{4};

%% compute empirical autocorrelation function 
m = size(y,1); R = zeros(taumax/dx,num);
for j = 2:num+1
    % Subtract mean from signal
    y(:,j) = y(:,j) - y(:,1);
    
    % R(i) = sum_k y_k*conj(y_{k-i})
    for i = 1:floor(taumax/dx)
        R(i,j) = y(i+1:m,j)'*y(1:m-i,j)/(m-i);
    end
end
R = mean(R,2)/sigma;