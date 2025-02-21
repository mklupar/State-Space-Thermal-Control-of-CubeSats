
clc;clear;close all
% Readme:

% Determines gains to produces desired system response.

% State Space parameters from SingelNodeThermalModel.m
A = [0.99967300061293151713925908552483];
B = [0.0010271976314499002731761567730473];
C = [1];
D = [0];
c = 0.073574862258581882890062786373164;
Ts = 1;

t1 = 1:Ts:30;
r1 = 300*ones(length(t1),1)';
r1(t1>=1) = 305;

% Augmented system to include constant
A_aug = [A ones(size(A,1));
    0 ones(size(A,1))];
B_aug = [B ; 0];
C_aug = [1 0; 0 1];
D_aug = 0;

sys_aug = ss(A_aug,B_aug,C_aug,D_aug,Ts);
pole = [0.999];
K1 = place(A,B,pole);
K2 = B\1;

G = inv(C*inv(1-(A-B*K1))*B); % Steady state tracking gain

A_cl = [(A-B*K1) eye(size(A)) - B*K2;
    0 eye(size(A))];
B_cl = [B*G ; 0];
C_cl = C_aug;
D_cl = D_aug;

sys_cl = ss(A_cl,B_cl,C_cl,D_cl,Ts);
t2 = 1:Ts:10000;
r2 = 300*ones(length(t2),1);
r2(t2>5) = 305;


figure(1)
[y,t]=lsim(sys_cl,r2,t2,[300 c]);
plot(t,y(:,1),'color','#0000FF')
title('Step Response')
ylabel('Temperature (K)')
xlabel('Time (s)')
grid on

h = figure(1)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'Step Response','-dpdf','-r0')

% Control Input
u = (G*r2 - K1*y(:,1) - K2*y(:,2));
figure(2)
plot(t,u,'color','#0000FF')
title('Input vs Time')
xlabel('Time (s)')
ylabel('Input (W)')
legend('Max Input = ' +string(max(u)) + 'W')
axis tight
grid on

h = figure(2)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'Step Input','-dpdf','-r0')
