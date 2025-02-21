
clc; clear; close all

% README: 
% kronDelta.m and TD_1Node.xls are needed for code to run. Files can be found in single-Node dir.
%(Note that TD_1Ndoe.xls will produces that same results if eight nodes were simulated given that they all have the same initial conditions)

% Define nodes and surfaces for Enclosure
numSurfaces = 9;
numNodes = 8;

% Dim of each control volume
dx = 0.1/2;
dy = 0.133/2;
dz = 0.3/2; 

% Surface area of each node
A = dx*dz + dy*dz + dx*dy;
A = A*ones(1,8);
Aenclosure = 0.4158;
A = [A Aenclosure];

%% Properties
epsilon = 0.5*ones(1,9);
T_inf = 2.73; % K
k_AL = 237; % W/m/K
rho_AL = 2702; % Denstiy of Aluminium (kg/m^3)
cp_AL = 903; % Specific heat of Aluminium (J/kg/K)
V = dx*dy*dz; % Volume 
alpha = 1/(rho_AL*cp_AL*V);
sigma = 5.6704E-8;

% View Factor
F = zeros(9,9);
F(1:end-1,end) = 1;

for i = 1:8
    F(9,i) = A(i)/Aenclosure;
end
F(end,end) = 1 - sum(F(end,:));

% Create nodes and identify surrounding node
n1 = [dx/2,dy/2,dz/2];
n2 = n1 + [0,0,dz];
n3 = n1 + [dx 0 0];
n4 = n1 + [dx 0 dz];
n5 = n1 + [0 dy 0];
n6 = n1 + [0 dy dz];
n7 = n1 + [dx dy 0];
n8 = n1 + [dx dy dz];

nodes = [n1; n2; n3; n4; n5; n6; n7; n8];

% Visulaize node locations
%figure
%scatter3(nodes(:,1), nodes(:,2), nodes(:,3))
%axis equal

surrNodeIndex = findSurrNodes(nodes,dx,dy,dz); % Row 1 is the nodes surrounding node 1, if there is no node in a given direction node 1 is placed at that index.

% Calculate radiation exchange (Enclosure Theory)
for i = 1:numSurfaces
    for j = 1:numSurfaces
        fluxCoeff(i,j) = (kronDelta(i,j)/epsilon(j) - (F(i,j)*(1-epsilon(j))/epsilon(j)))/A(j); % Flux Coefficients
        condCoeff(i,j) = (kronDelta(i,j)-F(i,j))*sigma; % Temperature Coefficients
    end
end

% [Q1 Q2 ... Q9 T1 T2 ... T9] (Interior Radiation)
M = [fluxCoeff -condCoeff];

%Shell Eq (External Radiation)
shellEq = zeros(1,numSurfaces);
shellEq(9) = 1;
shellEq(18) = epsilon(9)*sigma*A(9);

M = [M;shellEq];
M(:,end+1) = zeros(10,1);
M(end,end) = epsilon(9)*sigma*A(9)*T_inf^4;

% Reorganize and solve for T1-8 in terms of Q1_8 using rref
Mshell = [M(:,9) M(:,18)];
M(:,18) = []; M(:,9) = [];
%[Q9 T9 Q1 Q2 ..Q8 T1 T2 ... T8]
M = [Mshell M];
M = rref(M);

% Solve for temperaute of each node

T = 300*ones(numNodes,1);% initial temperature of each node
dt = 1;
endtime = 36000;
time1 = 0:dt:endtime;
u = 0*ones(8,endtime+1); % Internal heat generations (W)
% Nonlinear Calculations
for time = 0:dt:endtime
    
    tIndex = round(time/dt + 1);
    
    Qnonlinear(:,tIndex) = M(3:end,end) - M(3:end,11:end-1)*T(:,tIndex).^4;
   
    for i = 1:8
        node = surrNodeIndex(i,:); % Contains nodes surroudning node i

        T(i,tIndex+1) = alpha*k_AL*dt * ( dy*dz * ( T(node(1),tIndex) + T(node(2),tIndex) - 2*T(i,tIndex) )/dx ...
            + dx*dz * ( T(node(3),tIndex) + T(node(4),tIndex) - 2*T(i,tIndex) )/dy ...
            + dx*dy * ( T(node(5),tIndex) + T(node(6),tIndex) - 2*T(i,tIndex) )/dz) ...
            + T(i,tIndex) - Qnonlinear(i,tIndex)*alpha*dt + u(i,tIndex)*alpha*dt;
    end
end

TD_Results = table2array(readtable('TD_1Node.xls'));
time2 = TD_Results(:,1);
temperature = TD_Results(:,2);

figure(1)
plot(time2,temperature(:,:),'r')
hold on
plot(time1,T(:,1:end-1)','--','color','#0000FF')
hold off
grid on
axis tight
title ('System Enclosure - Nonlinear Radiation')
ylabel('Temperatuer (K)')
xlabel('Time (s)')

legend(["Thermal Desktop","MATLAB (Nonlinear)"])


h = figure(1)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'8Node Nonlinear','-dpdf','-r0')

% Linearize Radiation
Tsym = sym("T",[numNodes 1]);
Qnonlinear = M(3:end,end) - M(3:end,11:end-1)*Tsym.^4;
Qtaylor  = vpa(taylor(Qnonlinear,Tsym,300,Order=2)); % second order taylor series approx of radiation

% Convert equation into a matrix equation: T(t+1) = [A]T(t) + [B]u(t)
for row = 1:length(Qtaylor)
    radCoeff(row,:) = flip(coeffs(Qtaylor(row,:)));
end

Tlin = 300*ones(numNodes,1);
radCoeff = double(radCoeff);
const = radCoeff(:,end)*alpha*dt;
condCoeff = zeros(numNodes,numNodes);

for i = 1:numNodes
    xCoeff = alpha*k_AL*dt*dy*dz/dx;
    yCoeff = alpha*k_AL*dt*dx*dz/dy;
    zCoeff = alpha*k_AL*dt*dy*dx/dz;

    xNodes = surrNodeIndex(i,1:2);
    yNodes = surrNodeIndex(i,3:4);
    zNodes = surrNodeIndex(i,5:6);
    
    condCoeff(i,xNodes) = condCoeff(i,xNodes) + xCoeff;
    condCoeff(i,yNodes) = condCoeff(i,yNodes) + yCoeff;
    condCoeff(i,zNodes) = condCoeff(i,zNodes) + zCoeff;
end

for i = 1:numNodes
   condCoeff(i,i) = condCoeff(i,i) + (-2*alpha*k_AL*dt * ((dy*dz/dx) + (dx*dz/dy) + (dx*dy/dz)) + 1);
end

Taug = [T(:,1) ; -const(1)];

ssA = (condCoeff - alpha*dt*radCoeff(1:8,1:8));
ssB = alpha*dt*eye(numNodes);


ssA_aug = [ssA ones(8,1); zeros(1,8) 1];
ssB_aug = [ssB; zeros(1,8)];

for time = 0:dt:36000
    
    tIndex = round(time/dt + 1);

    Tlin(:,tIndex+1) = ssA*Tlin(:,tIndex) + ssB*u(:,tIndex) + (-const); % check to make sure no error in forming matrix equations
    Taug(:,tIndex+1) = ssA_aug*Taug(:,tIndex) + ssB_aug*u(:,tIndex);
    
     
end


figure(2)
plot(time2,temperature(:,1),'r')
hold on
plot(time1,Taug(1:8,1:end-1),"--",'color','#0000FF')
hold off
grid on
axis tight
title ('System Enclosure - Linear Radiation')
ylabel('Temperatuer (K)')
xlabel('Time (s)')
legend(["Thermal Desktop",'MATLAB (Linear)'])

h = figure(2)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'8Node Linear','-dpdf','-r0')

% save data to be imported into thermal controls file
save ssData.mat ssA ssB u dt const
