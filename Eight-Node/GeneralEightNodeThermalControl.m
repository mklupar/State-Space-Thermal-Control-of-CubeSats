% Thermal Control for a 8 node (m input, l output) model
clc; clear; close all
load ssData.mat

A = ssA;
B = ssB;

n = 8; % number of states
m = 8; % number of inputs: m can be changed to see the effects on steady state tracking as pinv() will be needed to account for V not being square.
l = 8; % number of outputs: In the case where the full state is not available at the output (L < 8), an observer will be needed (Kalman or Luenberger Observer).

B = B(1:n,1:m);

C = eye(l,n);
D = 0;
const = -const(1);
Ts = dt;

A_aug = [A ones(n,1);
    zeros(1,n) 1];
B_aug = [B; 
    zeros(1,m)];

C_aug = eye(l,n+1);
sys_aug = ss(A_aug,B_aug,C_aug,0,Ts);

pole = linspace(0.95,0.99,n);
K1 = place(A,B,pole);
K2 = B\ones(n,1);

I = eye(size(A));
V = C*inv(I-(A-B*K1))*B;
G = pinv(V);

A_cl = [A-B*K1 ones(n,1)-B*K2;
    zeros(1,n) 1];
B_cl = [B*G;
    zeros(1,l)];
C_cl = C_aug;
sys_cl = ss(A_cl,B_cl,C_cl,0,Ts)

t = 1:Ts:1000;
r = 300*ones(length(t),size(B_cl,2));
r(t>20,:) = 305;
init = [300*ones(1,n) const];

[y2,t2] = lsim(sys_cl,r,t,init);
figure(1)
plot(t2,y2)
title('Closed Loop System Response (7I/8O)')
xlabel('Time (s)')
ylabel('Temperature (K)')

% find states
state = [init(1:n)];
for i = t
    state(i+1,:) = ((A-B*K1)*state(i,1:n)' + (ones(n,1)-B*K2)*init(1,n+1:end)' + B*G*r(i,:)')'; % closed loop state transition 
end
state

h = figure(1);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'7I8O Node Closed Loop Step Response','-dpdf','-r0')



