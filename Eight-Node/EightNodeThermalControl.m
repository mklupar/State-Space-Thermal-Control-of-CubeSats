
clc; clear; close all

% README:
% Thermal Control for a 8 node (8 input, 8 output) model | B is square
% Data from EightNodeThermalModeling.m

% (A, B) from Finite Difference Equations
load ssData.mat
A = ssA;
B = ssB;

C = eye(8,8);
D = 0;
const = -const(1);
Ts = dt;

% Augment system by adding constant term as an uncontrolled state
A_aug = [A ones(8,1);
    0 0 0 0 0 0 0 0 1];
B_aug = [B; 
    zeros(1,8)];
  
C_aug = eye(8+1,8+1);

sys_aug = ss(A_aug,B_aug,C_aug,0,Ts) % form state space representation

% Determine gains from desired pole placements
pole = [0.8 0.83 0.86 0.9 0.93 0.96 0.98 0.99] ;
Ac = sys_aug.A(1:8,1:8);
A12 = sys_aug.A(1:8,9);
Auc = sys_aug.A(9,9);

k1 = place(Ac,B,pole); % Gain for controllable subsystem
k2 = B\A12; % gain for uncontrollable subsystem
G = inv(C*inv(eye(8)-(A-B*k1))*B); % Steady state tracking gain

% Form closed loop system
Anew = [(Ac - B*k1) A12-B*k2;
    zeros(1,8) Auc];

Bnew = [B*G; zeros(1,8)];

Cnew = eye(9);
sys_new = ss(Anew,Bnew,Cnew,[],Ts)

t = 1:Ts:500;
r = 305*ones(length(t),8);
r(t<5,:) = 300;
init = [300*ones(1,8) const];
[y2,t2] = lsim(sys_new,r,t,init);

% Plot results
figure(1)
plot(t2,y2(:,1:8))
title('Closed Loop System Response')
xlabel('Time (s)')
ylabel('Temperature (K)')

h = figure(1)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'8 Node Closed Loop Step Response','-dpdf','-r0')

% required input

u = (G*r' - k1*y2(:,1:8)' - k2*y2(:,end)');

figure(2)
plot(t2,u')

h = figure(2)

title('Input Response')
ylabel('Input (W)')
xlabel('Time (s)')
legend('Max Input = ' +string(max(u(1,:))) + ' W')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'8 Node Closed Loop Input Response','-dpdf','-r0')


