function procAB_samplePaths(par)

%% parameters
T = par{1};
L = par{2};
v0 = par{3};
gamma = par{4};
N = par{5};
folder = par{6};
deter = par{7};
saveInt = par{8};
numTraj = 100;
l0 = v0; %average run length
x0 = L/3;

%% load data
disp(gamma)
if deter == 1; ind = 0; elseif deter == 0; ind = 1:Nprof; end
file = ['/AB_T_',num2str(T),'_N_',num2str(N),...
    '_gamma_',num2str(gamma),'_',num2str(ind),'.mat'];
if exist([folder file], 'file') == 0 
    disp([folder file])
    error('No result file!')
else
    data = load([folder file]);
    X = data.X;         
end

%% plot sample paths
figure; hold on
y = linspace(1,T,T/saveInt+1);
for i = 1:numTraj
    x = (X(i,:) - x0)/l0;
    z = zeros(size(x));
    surface([x;x],[y;y],[z;z],'facecol','no','edgecol','interp',...
                              'linew',0.5,'edgealpha',.1,'edgecolor','b');  
end

box(gca,'off')
set(gca,'FontSize',20)
xlabel('x [l_0]','FontSize',25)
set(gca,'XTick',[-2 0 2 4 6]/l0)
xlim([-2 4]/l0)
ylabel('time [\lambda_0^{-1}]','FontSize',25)
ylim([0 T])
set(gca,'YTick',linspace(0,T,3))

%% mean trajectories

%empirical drift
meanEmp = -x0/l0 + mean(X,1)/l0;

%exact
D = v0^2/2;
B = data.beta;
vexact = 2*D*B.*gamma*gamma./(1 + 2*gamma).^3;
meanExact = vexact.*y;

%% plot
plot(meanEmp,y,'LineWidth',2,'Color',[0 0 0],'LineStyle','-');
plot(meanExact/l0,y,'LineWidth',2,'Color',[1 0 0],'LineStyle','--');