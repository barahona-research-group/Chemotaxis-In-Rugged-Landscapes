function testKt
% This script tests if the impulse response of the dynamical system gives 
% consistent results with the tumbling kernel K(t).
%
% conclusion : for gamma>0.25, need dt <=0.005
%              for 0.1<gamma<0.25, need dt <=0.001                

%% parameters
gamma = 1;
beta = 1;
dt = 0.001;
T = gamma*10; %simulate for 10 memory lengths

%delta input
St = ones(1,floor(T/dt));
St(2) = 1/dt;

%initial conditions
m(1) = St(1)*gamma;
m(2) = St(1)*gamma^2;
m(3) = 2*St(1)*gamma^3;

%compute response
for i = 1:T/dt
   [m, Lambda(i)] = Kt(m,St(i),gamma,dt);
end

% exact Kt
fac = exp(sqrt(2)-2)*(sqrt(2)-1);
K = @(t) beta/fac*exp(-t/gamma).*(t/gamma - t.^2/(2*gamma^2));

% plot comparison
t = dt*(dt:floor(T/dt));
plot(t,Lambda/fac)
hold on
plot(t,K(t))