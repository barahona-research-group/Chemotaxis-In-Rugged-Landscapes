function [m, Lambda] = Kt(m,S,gamma,dt)

%This script implements the linear dynamical system modelling the
%chemotactic response K(t) using a first-order Euler scheme

%internal dynamics
m(:,1) = m(:,1) + dt * (        S - m(:,1)/gamma );
m(:,2) = m(:,2) + dt * (   m(:,1) - m(:,2)/gamma );
m(:,3) = m(:,3) + dt * ( 2*m(:,2) - m(:,3)/gamma );
    
%tumbling response
Lambda = m(:,2)/gamma - m(:,3)/(2*gamma^2);