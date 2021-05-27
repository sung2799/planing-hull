% Single run simulations, with mass & rudder options 3-30-21

%Initial Condition: 
IC=1; %1-Straight trajectory, 2-SS of previous simulation
if IC==1
clearvars -except IC; 
% close all;  % 4DOF, SISO Single run, Gain Sched interp model, measure axes, custom integrator 1/6/21
end

%alerts
alert=2; %alert: 1-on, 2-off
SD=2; %system shutdown after sim: 1-on, 2-off

%saving
SF=2; %save file: 1-on, 2-off
SG=2; %save graphs: 1-on, 2-off

%Simulation parameters
int=2; %Interpolation method: 1-Innate (Slow), 2-Expedient (Fast)
tend=200; 

%actuator inputs
I=2; %Actuator input: 1-mass, 2-rudder or both
U_m=-.3; d0=U_m; %U_mm=0.3; 
U_r=-3; %U_r=-3.035  DR10: reference rudder input to match yaw rate of U_m=.3

%Geometric parameters
lcg=1.2; %range, .8(t0=6) - 1.5(t0=3) DR10, .97(t0=6) - 1.63(t0=3) DR20
DR=20; % Deadrise

sim_setup

U_x=2;

[XX,YY,tri]= interp_model(lcg);

% % enter conditions to test & linearize
if IC==1
x0=[0;0;0;0;14;0;0;0];  %IC
end

if IC==2
x0=yNL(end,:)';  %steady state condition
end

%   x y F S u  v p r
x1=x0;  %allows setting of an initial perturbation or reference tracking

U0=[U_x;U_m;U_r];

% Make controller
% [AA,BB,KK,K2] = SISO_LQR(x0,XX,YY,tri,lcg,vcg,DR,ma,mb,A,B,t0,U0,I,int); 
% KK(:,4)=0; % remove S driver
% [V,E,W]=eig(AA); %Columns of W are left eigenvectors ST W'*A = D*W'
% CC=W'*BB;
% CC_m=norm(CC);

% Linearize to AA matrix system only
[AA,BB] = SISO_LQR(x0,XX,YY,tri,lcg,vcg,DR,ma,mb,A,B,t0,U0,I,int); 
KK(:,4)=0; % remove S driver
[V,E,W]=eig(AA); E=diag(E); %Columns of W are left eigenvectors ST W'*A = D*W'
PMM=AA(3:end,3:end)*A;
CC=W'*BB;
CC_m=norm(CC);


%% Simulation
%yNL: non-linear response
dt=.01; 
tspan = dt:dt:tend; iter1=size(tspan,2);
tstart=tic;

[t,yNL] = ode15s(@(t,y)roll_dynamics_4DOF_NL16(y,XX,YY,tri,lcg,vcg,lk,b,Cv,DR,ma,mb,A,B,t0,U0,int,I),tspan,x0); %Open Loop 

tr=toc(tstart);

%% Plots
figure(1) % Non linear Graph
title("Simulated Trajectories, "+ DR +"^{o} Deadrise ");
subplot(611), plot(t,yNL(:,3)); grid;
ylabel('\phi(^o)');
hold on
subplot(612), plot(t,yNL(:,4)); grid; 
ylabel('\psi(^o)');
hold on
subplot(613), plot(t,yNL(:,5)); grid; 
ylabel('u(ft/s)');
hold on
subplot(614), plot(t,yNL(:,6)); grid; 
ylabel('v(ft/s)');
hold on
subplot(615), plot(t,yNL(:,7)); grid; 
ylabel('p(^o/s)');
hold on
subplot(616), plot(t,yNL(:,8)); grid; 
ylabel('r(^o/s)');
hold on
xlabel('Time(s)')
% legend('baseline, No \phi','\tau=6^{o}, +\phi','\tau=3^{o}, +\phi')
% legend('Mass actuated','Rudder actuated')

figure(2) % x-y plot
title("Trajectories in x-y, "+ DR +"^{o} Deadrise")
plot(yNL(:,1),-yNL(:,2))
xlabel('ft')
ylabel('ft')
hold on
% legend('Mass actuated','Rudder actuated')

figure(3)
plot(t,yNL(:,8)); grid;
title('Yaw rate')
xlabel('Time(s)');
ylabel('r');
hold on
% legend('Mass actuated','Rudder actuated')

% Control mass position history if actively controlled
% U_m_his=KK*(yNL'-x0)+U0(1:ci);
% figure(4)
% plot(t,U_m_his)
% xlabel('time(s)')
% title('Control input history (MIMO)')
% legend('U_x','U_m','U_r')
% hold on

% Final power metrics
u=yNL(:,5); v=yNL(:,6); F=yNL(:,3); S=yNL(:,4); r=yNL(:,8);
uw_f=(u(end)^2+v(end)^2)^.5;
XI=[DR t0 v(end) r(end) F(end)];
[Yint]=interp_point(XX,YY,tri,XI,int); 
Xf=Yint(1); Yf=Yint(2); Kf=Yint(3); Nf=Yint(4);
Xw_f=(Xf^2+Yf^2)^.5;  %wind axes forces
Pw_f=uw_f*Xw_f; %Power history
AoA_f=U_r-atand(v(end)/u(end));
if I==1
    Xrud_f=0;
end
if I==2
    Xrud_f=-.5*(4.6e-5*AoA_f^2+1.38e-4*abs(AoA_f)+9.2e-4)*u(end)^2;
end
Pr_f=u(end)*(-Xrud_f);
Pt_f=Pr_f+Pw_f
r(end)

% Final turning metrics
rt=(u(end)^2+v(end)^2)^.5*360/(2*pi*r(end));
Ad=max(yNL(:,1));
TD=max(yNL(:,2));
AdTD=Ad/TD;

% %
if SG==1
saveas(figure(1),'sim123a');saveas(figure(1),'sim123b.jpg');
saveas(figure(2),'sim123c');saveas(figure(2),'sim123d.jpg');
saveas(figure(3),'sim123e');saveas(figure(3),'sim123f.jpg');
end

if SF==1
save('sim_data.mat');
end

% Video

%% Reproduction of forces
% Rudder modeling
if I==1
u=yNL(:,5); v=yNL(:,6); F=yNL(:,3); S=yNL(:,4); r=yNL(:,8);
C_X=0; %setting thrust = drag for constant velocity
C_Y=0; %sway from rudder
C_K=0; % roll force from control mass & steering mechanism
C_N=0;
end

if I==2
u=yNL(:,5); v=yNL(:,6); F=yNL(:,3); S=yNL(:,4); r=yNL(:,8);
z_rud=vcg; %ft

%Simplified rudder model
SL=atand(v./u);
AoA=U_r-SL;
L=.5*(.0018)*u.^2.*AoA;
Y_rud=L;
X_rud=-.5*(4.6e-5*AoA.^2+1.38e-4*abs(AoA)+9.2e-4).*u.^2;
K_rud=-Y_rud.*z_rud;
N_rud=-Y_rud.*lcg;

% Group Control forces
C_X=U_x+X_rud; %setting thrust = drag for constant velocity
C_Y=Y_rud; %sway from rudder
C_K=ma*g*U_m+K_rud; % roll force from control mass & steering mechanism
C_N=N_rud; %thrust yaw coupling likely double counted in HD forces
end

%% Label & Save

if SG==1
saveas(figure(1),'sim123a');saveas(figure(1),'sim123b.jpg');
saveas(figure(2),'sim123c');saveas(figure(2),'sim123d.jpg');
saveas(figure(3),'sim123e');saveas(figure(3),'sim123f.jpg');
saveas(figure(5),'sim123g');saveas(figure(5),'sim123h.jpg');
end

if SF==1
save('sim_data.mat');
end

if alert==1
load handel;
sound(y)
end
     
if SD==1
% system('shutdown -s')
system('sleep -s')
end