% Lax-Friedrichs method to solve 1-D Euler equations Isothermal
% Modificado por Julio Cesar Santos Nascimento
% Data: 26/06/2014  08:00 as 11:00

clear;
n = 101;         %Number of grid points
L = 4;           %Length of domain
h = L/(n-1);      %Spatial step size
CFL = 0.34;      %CFL number for stability
t_final = 2e-3;  %Final time
x =0:h:L;        % x grid position
R_ar = 287.038;
T = 450;
%-------------------------------------------------------------------------%
% Define intial conditions 
%-------------------------------------------------------------------------%
p_l = 500e3;  %Pressure in left side of shock tube at t=0
p_r = 100e3;  %Pressure in right side of shock tube at t=0
rho_l = 3.87;   %Density at left side of shock tube at t=0
rho_r = 0.774; %Density at right side of shock tube at t=0
u_l = 0;     %Velocity in left side of shock tube at t=0
u_r = 0;     %Velocity in right side of shock tube at t=0

p(1:(n+1)/2) = p_l;
p((n+3)/2:1:n) = p_r;
Temp(1:(n+1)/2) = T;
Temp((n+3)/2:1:n) = T;
rho(1:(n+1)/2) = rho_l;
rho((n+3)/2:1:n) = rho_r;
u(1:(n+1)/2) = u_l;
u((n+3)/2:1:n) = u_r;
  %Total Energy
a = sqrt(rho*R_ar*T);            %Speed of sound
dt = CFL*h/max(abs(u+a));          %Time step based on CFL number
step = 0;
%-------------------------------------------------------------------------%
% Time integration begins 
%-------------------------------------------------------------------------%
for t = dt:dt:t_final
    % Define q & F matrix
    q = [rho; rho.*u];
    F = [rho.*u; rho.*(u.^2)+p];
    % Update q matrix and Flow Parameters
    q(1:2,2:n-1) = 0.5*(q(1:2,3:n)+q(1:2,1:n-2))-(dt/(2*h))*(F(1:2,3:n)-F(1:2,1:n-2));
    rho = q(1,1:n);
    u = q(2,1:n)./rho;
    p = rho.*R_ar*T;
    Temp = p/(rho.*R_ar);
    step = step+1;
end

% Calculation of Flow Parameters
a = sqrt(rho*R_ar*T); 
M = u./a;
p_ref = 101325;
rho_ref = 1.225;
Q = rho.*u;
%-------------------------------------------------------------------------%
% Plot the variables 
%-------------------------------------------------------------------------%
figure 
offset = 0.05;
subplot(2,3,1);plot(x,p,'-k'); xlabel('X-Coordinate (m)');ylabel('Pressure (Pa)');

subplot(2,3,2);plot(x,u,'-k'); xlabel('X-Coordinate (m)');ylabel('Velocity (m/s)');

subplot(2,3,3);plot(x,rho,'-k'); xlabel('X-Coordinate (m)');ylabel('Density (kg/m³)');

subplot(2,3,4);plot(x,M,'-k'); xlabel('X-Coordinate (m)');ylabel('Mach Number ');

subplot(2,3,5);plot(x,Q,'-k'); xlabel('X-Coordinate (m)');ylabel('Mass Flow (kg/s ');

subplot(2,3,6);plot(x,Temp,'-k'); xlabel('X-Coordinate (m)');ylabel('Temperature k) ');

%-------------------------------------------------------------------------%
