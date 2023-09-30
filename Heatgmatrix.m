% Solves the 1D heat equation with an explicit finite difference scheme
% 
% From Thorsten Becker's webpage
%

clear all
close all

% Physical parameters
L = 100; % Length of modeled domain [m]
Tmagma = 1200; % Temperature of magma [C]
Trock = 300; % Temperature of country rock [C]
kappa = 1e-6; % Thermal diffusivity of rock [m2/s]
W = 5; % Width of dike [m]
day = 3600*24; % # seconds per day
dt = 1.4*day; % Timestep [s]

% Numerical parameters
nx = 201; % Number of gridpoints in x-direction
nt = 100; % Number of timesteps to compute
dx = L/(nx-1); % Spacing of grid
x = -L/2:dx:L/2;% Grid

% Setup initial temperature profile
T = ones(size(x))*Trock;
T(find(abs(x)<=W/2)) = Tmagma;

time = 0;
c=(kappa*dt)/(dx*dx);
G=zeros(nx,nx);
 %TG=G*T

% Compute new temperature
    for i=1:nx
        if(i==1||i==nx)
            G(i,i)=1;
        else 
            G(i,i+1)=c;
            G(i,i)=1-2*c;
            G(i,i-1)=c;
        end
    end
    T2=transpose(T);
  for n=1:nt % Timestep loop
       Tnew=G*T2;
       T2=Tnew;
       time=time+dt;
    
    % Plot solution
    figure(1), clf
    plot(x,Tnew);
    Ymat(n,:) = n*dt*(ones(1,nx));
    Xmat(n,:) = x;
    Tmat(n,:) = Tnew;
    axis([-50 50 200 1300])
    xlabel('x [m]')
    ylabel('Temperature [^oC]')
    title(['Temperature evolution after ',num2str(time/day),' days'])
    drawnow
end
figure
surf(Xmat,Ymat,Tmat)
