
clc;clear;close all; warning 'off'

% README:

% Requires kronDelta.m and TD_1Node.xls to be in cwd.

% Lumped system model for CubeSat. Optics box is modeled as a solid
% box (represented as a node) enclouse by the structure (represented as
% a thin shell box).



%% Properties

% Optics Box dim
dx = 0.1; %m
dy = 0.3; %m
dz = 0.133; %m

% Structure dim
DX = 0.2263; %m
DY = 0.2263; %m
DZ = 0.346; %m

sigma = 5.6704E-8; % Stefan Boltzman Const. (W/m^2/K^4)

A1 = 2*dx*dy + 2*dx*dz + 2*dy*dz; % Surface area of optics box (m^2)
A2 = 2*DX*DY + 2*DX*DZ + 2*DY*DZ; % Surface area of enclosing structure (m^2)
Vs = dx*dy*dz; % Volume of optics box (m^3)
A = [A1 A2];

eps1 = 0.5; % Emisivity of optics box
eps2 = 0.5; % Emissivity of enclosing structure
eps = [eps1 eps2 eps2];

F = [0 1 ; A1/A2 1-A1/A2]; % View Factor matrix

T_inf = 2.73; % Ambeint Temperature (K)

rho_AL = 2702; % Denstiy of Aluminium (kg/m^3)
cp_AL = 903; % Specific heat of Aluminium (J/kg/K)

%% Surface matrix equations (Enclosure Theory)
for i = 1:2
    for j = 1:2
        E(i,j) = (kronDelta(i,j)/eps(j) - (F(i,j)*(1-eps(j))/eps(j)))/A(j); % Flux Coefficients
        D(i,j) = (kronDelta(i,j)-F(i,j))*sigma; % Temperature Coefficients
    end
end

% Solve Enclosure systems of equations for Q1 as that is the net radiation loss term
% [Q1 Q2 T1 T2]
M = [E -D];
M(3,:) = [0 1 0 eps2*sigma*A2]; % Add 3rd equation for external radiation
M(:,5) = [0 0 eps2*sigma*A2*T_inf^4]'; 
%Reorgamize matrix to be [T2 Q2 T1 Q1]
M = [M(:,4) M(:,2) M(:,3) M(:,1) M(:,5)]; 
M = rref(M); % Use rref() to 'simlify/solve' system of equations

%% Nonlinear transient solution for T1
dt = 1; % time step (s)
endtime = 36000; % time of experiment (s)
time_index = 1:endtime/dt;
time = 0:dt:endtime;

T1 = [300]; % initial temperature of node 1 (K)
Q1 = [];
u = 0; % input into optics (W)

for t = time_index
    Q1(end+1) = -(T1(t)^4 - M(3,5))/M(3,4);
    T1(end+1) = -Q1(t)*dt/(rho_AL*cp_AL*Vs) + u*dt/(rho_AL*cp_AL*Vs) + T1(t); % Solve "E_st = E_in - E_out + E_gen" for T1(t+1)
    
end

%% Nonlinear Plots
% Load TD results
TD_Results1 = table2array(readtable('TD_1Node.xls'));
time1 = TD_Results1(:,1);
temperature1 = TD_Results1(:,2);

figure(1)
plot(time1,temperature1, 'r')
hold on
plot(time,T1',"--",'color','#0000FF')
hold off
title('System Enclosure - Nonlinear Radiation')
legend('Thermal Desktop','MATLAB Nonlinear')

xlabel('Time(s)')
ylabel('Temperature (K)')
grid on 
axis tight

h = figure(1)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'System Enclosure - Nonlinear Radiation','-dpdf','-r0')
%% Linear trasnient solution for T1

T1_Linear = [300]; % temp of node 1 (K)
T_NoShift = [300];
Q1_Linear = [];
operating_point = 300;
for t = time_index
    
    % 2 term taylor series approx of nonlinear Q1 around operating point
    %Q1_Linear(end+1) = ((-operating_point^4 + M(3,5))/M(3,4) - (4*operating_point^3 / M(3,4))* (T1(t)-operating_point));
    % Linear Temperature equation
    %T1(end+1) = -Q1(t)*dt/(rho_AL*cp_AL*Vs) + u*dt/(rho_AL*cp_AL*Vs) + T1(t);
    
    % Sub Q1 into T1
    T_NoShift(end+1) = (((4*operating_point^3*dt)/(M(3,4)*rho_AL*cp_AL*Vs)) + 1)*T_NoShift(t) +  (dt/(rho_AL*cp_AL*Vs))*u ;%-((3*operating_point^4 + M(3,5))*dt)/(M(3,4)*rho_AL*cp_AL*Vs);
    T1_Linear(end+1) = (((4*operating_point^3*dt)/(M(3,4)*rho_AL*cp_AL*Vs)) + 1)*T1_Linear(t) +  (dt/(rho_AL*cp_AL*Vs))*u -((3*operating_point^4 + M(3,5))*dt)/(M(3,4)*rho_AL*cp_AL*Vs);
end

% Exam shift by const.
% Control can be done with out without the sift from the const. term, this
% just shows that everything matches.
e = [0];
A_m = (((4*operating_point^3*dt)/(M(3,4)*rho_AL*cp_AL*Vs)) + 1);
shift = -((3*operating_point^4 + M(3,5))*dt)/(M(3,4)*rho_AL*cp_AL*Vs);
for i = 1:t
    e(end+1) = e(i) + A_m^i * shift;
end

T_Shifted = T_NoShift+e;

%% Linear Plots 
figure(2)
plot(time1,temperature1,'color','r')
hold on
plot(time,T1_Linear,'--','color','#0000FF')

title('System Enclosure - Linear Radiation')
legend('Thermal Desktop','MATLAB Linear')
grid on 
axis tight
%axis([0 2000 295 300])

xlabel('Time (s)')
ylabel('Temperature (K)')

h = figure(2)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'System Enclosure - Linear Radiation','-dpdf','-r0')

figure()
plot(time,T1_Linear,'color','r')
hold on
plot(time,T_NoShift,'color','b')
plot(time,T_Shifted,'--','color',"k")
hold off
title('Shifted Dynamics')
xlabel('Time')
ylabel('Temperature (K)')
legend('Ax + Bu + \phi','Ax + Bu', 'Ax + Bu + shift')

%% Steady State Results (Just simulating for a longer period to show the effects of linearization)

dt = 10; % time step (s)
endtime = 3600000; % time of experiment (s) | just make time go big
time_index = 1:endtime/dt;
time = 0:dt:endtime;

T1 = [300]; % temp of node 1 (K)
T1_Linear = [300];
Q1 = [];
u = 0; % input into optics (W/m^2)
for t = time_index
    Q1(end+1) = -(T1(t)^4 - M(3,5))/M(3,4);
    T1(end+1) = -Q1(t)*dt/(rho_AL*cp_AL*Vs) + u*dt/(rho_AL*cp_AL*Vs) + T1(t);  
    T1_Linear(end+1) = (((4*operating_point^3*dt)/(M(3,4)*rho_AL*cp_AL*Vs)) + 1)*T1_Linear(t) +  (dt/(rho_AL*cp_AL*Vs))*u -((3*operating_point^4 + M(3,5))*dt)/(M(3,4)*rho_AL*cp_AL*Vs);
end
%%Steady State Plots
figure()
plot(time,T1,'r')
hold on
plot(time,T1_Linear,'color','#0000FF')
hold off

title('Steady State Results')
xlabel('Time (s)')
ylabel('Temperature (K)')
legend('Nonlinear', 'Linear')
grid on 
axis tight
