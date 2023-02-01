clc;clear;close all;
%% Given Properties
Q = 0.6; % Heat Generation Rate [cal/cm^3s]
k = 0.16; % Thermal Conductivity [cal/cm°C]
H = 0.073; % Heat Transfer Coefficient [cal/cm^2°C]
tau = 0.5; % Plate Thickness [cm]
T0 = 25; % Surrounding Temperature [°C]
%% Preallocation
xlength = 9; % x length [cm]
ylength = 5; % y length [cm]
h = 0.25; % cm, Interval
xnode = xlength/h + 1; % Number of Nodes on x axis
ynode = ylength/h + 3; % Number of Nodes on y axis
T = zeros(ynode,xnode); % Final Preallocated Array for Temperature
%% SOR Factor Calculations
syms w
p = xlength/h; % x size divided by step size
q = ylength/h; % y size divided by step size
F = (((cos(pi/q))+(cos(pi/p)))^2)*w^2-16*w+16==0; % SOR factor equation
w = double(vpasolve(F,w)); % Solution of the equation
for i=1:1:length(w)
    if 0<w(i) && w(i)<2 % SOR factor changes between 0 and 2
        omega=w(i);
    end
end
%% Drichlet Boundary Conditions
Tleft = 20; % Left Side Temperatures
Tright = Tleft; % Right Side Temperatures
T(2:end-1,1) = Tleft;
T(2:end-1,end) = Tright;
%% SOR Iteration Calculations
for kiter = 1:100000000
    a = T(2,2);
    for i=2:xnode-1
        % Lower Edge Points
        T(end,i) = T (end-2,i)-30*h; % Lower Imaginary Line
        T(end-1,i) = T(end-1,i) + omega*((T(end-1,i-1)+T(end,i)+T(end-1,i+1)+T(end-2,i)-4*T(end-1,i))/4 + (Q*h^2)/(4*k*tau)); % Lower Boundary Line
        % Upper Edge Points
        T(1,i) = T(3,i) - 2*h*H/k*(T(2,i)-T0); % Upper Imaginary Line
        T(2,i) = T(2,i) + omega*((T(2,i-1)+T(3,i)+T(2,i+1)+T(1,i)-4*T(2,i))/4 +(Q*h^2)/(4*k*tau)); % Upper Boundary Line
        % Interior Points
        for j = 2:ynode-1
            T(j,i) = T(j,i) + omega*((T(j,i-1)+T(j-1,i)+T(j,i+1)+T(j+1,i)-4*T(j,i))/4 + (Q*h^2)/(4*k*tau));
        end
    end
    b = T(2,2);
    err = abs((a-b)/b)*100;
    if err < 10^-6
        break
    end
end
%% Heat Map Plotting
% Temperature Distribution with Imaginary Nodes
figure(1)
Timaginary = T;
heatmap(T,'colormap',turbo);
title('Temperature Distribution with Imaginary Nodes');
xlabel('x = 9 cm');
ylabel('y = 5 cm');
% Temperature Distribution without Imaginary Nodes
figure(2)
T = T(2:end-1,:);
heatmap(T,'colormap',turbo);
title('Temperature Distribution without Imaginary Nodes');
xlabel('x = 9 cm');
ylabel('y = 5 cm');