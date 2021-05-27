% Rudder & Mass actuator on LQR controller on 5DOF, new lagrangian dynamics, 5-7
clear all;close all;clc
load DR20_lcg15_Ur4_Um-15_MIMO_514  %new actuator dynamics, mass fixed toggle so AA is stable 5-14

% Provides linear dynamics (AA,BB) about any given state "x0" with ENV disturbance channel in BB
% [AA,BB] = MIMO(x0,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,I,int,CL);

%saving
SF=0; %save file: 1-on, 0-off
SG=0; %save graphs: 1-on, 0-off

%LQR Controller development
Q=diag([0.4789;  %using least stable oscillatory OL eigenvector for 5DOF system
    0.0319;
    0.0112;
    0.0227;
    0.8442;
    0.0562;
    0.1138;
    0.2006]);

Q0=diag([-0.2938;  %using least stable oscillatory OL eigenvector for 4DOF system
    0.9333;
   -0.0570;
    0.0032;
    0.0597;
   -0.1895]);

% Make LQR controllers & CL eigs
R2=[1 0; 0 1]; R1=1;
AA0=AA0(1:6,1:6);  %mass fixed 4DOF system
B1=BB(:,1); B2=BB(:,1:2); B1r=BB(1:6,2);
K1=lqr(AA,B1,Q,R1);
K2=lqr(AA,B2,Q,R2);
K1r=lqr(AA0,B1r,Q0,R1);
[T0,E0]=eig(AA0); E0=diag(E0);
[T1,E1]=eig(AA-B1*K1); E1=diag(E1);
[T2,E2]=eig(AA-B2*K2); E2=diag(E2);
[T1r,E1r]=eig(AA0-B1r*K1r); E1r=diag(E1r);
B_di=BB(:,3); %create disturbance input channel on state F
B_di0=BB(1:6,3); %create disturbance input channel on state F
% C=[AA(5,:);0 0 0 0 0 0 1 0]; D=zeros(1,1); %output: pdot, d1

%Pole Zeros
C2=[1 0 0 0 0 0 0 0;0 0 0 0 0 0 1 0]; D=zeros(1,1); %output: F, d1
C0=[1 0 0 0 0 0]; %F
OL_sys=ss(AA0,B_di0,C0,D);  %open loop LQR (mass fixed) with disturbance channel
LQR_sys1=ss(AA-B1*K1,B_di,C2,D); %closed loop single actuator w/ disturbance channel
LQR_sys2=ss(AA-B2*K2,B_di,C2,D); %closed loop double actuator w/ disturbance channel
LQR_sys1r=ss(AA0-B1r*K1r,B_di0,C0,D); %closed loop double actuator w/ disturbance channel
figure(6)
pzplot(OL_sys); hold on
pzplot(LQR_sys1); hold on
pzplot(LQR_sys2); hold on
pzplot(LQR_sys1r); hold on
legend('Open Loop - mass fixed','Single Actuator, AMA','Double Actuator, AMA + RRS','Single Actuator, RSS & mass fixed')
% saveas(figure(1),'PZ');saveas(figure(1),'PZ.jpg');

%% Linear state space systems
% C2=AA(5,:);  %measure pdot
K1=[K1(:,1) zeros(1,2) K1(:,4:end)];
K2=[K2(:,1) zeros(2,2) K2(:,4:end)];
K1r=[K1r(:,1) zeros(1,2) K1r(:,4:end)];
C2=[1 0 0 0 0 0 0 0]; 
%   F S u v p r d1d2
sys0=ss(AA0,B_di0,C0,D);
sys1=ss(AA-B1*K1,B_di,C2,D);
sys2=ss(AA-B2*K2,B_di,C2,D);
sys1r=ss(AA0-B1r*K1r,B_di0,C0,D);

%Bode
figure(7)
bode(sys0); hold on;
bode(sys1); hold on;
bode(sys2); hold on;
bode(sys1r);
legend('Open Loop - mass fixed','Single Actuator, AMA','Double Actuator, AMA + RRS','Single Actuator, RSS & mass fixed')
figure(8)
step(sys1); hold on;
step(sys2); hold on;
step(sys1r);
title('Linear response to step input')
legend('Open Loop - mass fixed','Single Actuator, AMA','Double Actuator, AMA + RRS','Single Actuator, RSS & mass fixed')

% saveas(figure(2),'bodeMIMO');saveas(figure(2),'bodeMIMO.jpg');

%% Time domain linear sim
dt=.0025;
t=0:dt:3;
dist=zeros(size(t));
dist(1:end)=dt*(sin(8*pi*t(1:end)));  %generate wave disturbance, 1:101 if single impulse
y0=lsim(sys0,dist,t);  %open loop system response to wave disturance
y1=lsim(sys1,dist,t);  %closed loop system, 1 actuator response to wave disturance
y2=lsim(sys2,dist,t);  %closed loop system, 2 actuator response to wave disturance 
y1r=lsim(sys1r,dist,t);
figure(9)
plot(t,y1,t,y2,t,y1r,t,dist,'y')
legend('Single Actuator, AMA','Double Actuator, AMA + RRS','Single Actuator, RSS','Disturbance')
ylim([-.003 .003])
xlabel('Time(s)')

%% NL Sim

% Open loop
CL=0;
tend=100;
dt=.01; 
tspan = dt:dt:tend; 
x1=x0+0.01*[zeros(2,1);abs(T0(:,4));zeros(2,1)];
[t,yNL] = ode15s(@(t,y)roll_dynamics_4DOF_NL18(y,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,int,I,CL),tspan,x1); %Closed Loop
% [t,yNL] = ode15s(@(t,y)roll_dynamics_4DOF_NL18(y,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,int,I,CL),tspan,x0); %Open Loop 

Draw

% Single actuator
KK=[zeros(1,2) K1]; %LQR, fit 8x8 system to 10x10 (LQR)
KK(:,4)=0; %zero out control for heading
KK(:,5)=0; %zero out control for u velocity
KK(:,6)=0; %zero out control for v velocity
x1=x0+0.01*[zeros(2,1);abs(T1(:,4))];
CL=1;
[t,yNL] = ode15s(@(t,y)roll_dynamics_4DOF_NL18(y,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,[U0(1);KK*(x0-y);U0(3:4)],int,I,CL),tspan,x1); %Closed Loop
% [t,yNL] = ode15s(@(t,y)roll_dynamics_4DOF_NL18(y,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,int,I,CL),tspan,x0); %Open Loop 

Draw

% Double actuator
KK=[zeros(2,2) K2]; %LQR, fit 8x8 system to 10x10 (LQR)
KK(:,4)=zeros(2,1); %zero out control for heading
KK(:,5)=zeros(2,1); %zero out control for u velocity
KK(:,6)=zeros(2,1); %zero out control for v velocity
x1=x0+0.01*[zeros(2,1);abs(T2(:,4))];
CL=1;
[t,yNL] = ode15s(@(t,y)roll_dynamics_4DOF_NL18(y,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,[U0(1);KK*(x0-y);U0(4)],int,I,CL),tspan,x1); %Closed Loop  
% [t,yNL] = ode15s(@(t,y)roll_dynamics_4DOF_NL18(y,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,int,I,CL),tspan,x0); %Open Loop 

Draw

% RRS
KK=[zeros(1,2) K1r zeros(1,2)]; %LQR, fit 8x8 system to 10x10 (LQR)
KK(:,4)=0; %zero out control for heading
KK(:,5)=0; %zero out control for u velocity
KK(:,6)=0; %zero out control for v velocity
x1=x0+0.01*[zeros(2,1);abs(T0(:,4));zeros(2,1)];
CL=0; %Freeze position of control mass = mimics RRS with asymmetrical mass
[t,yNL] = ode15s(@(t,y)roll_dynamics_4DOF_NL18(y,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,[U0(1:2);KK*(x0-y);U0(4)],int,I,CL),tspan,x1); %Closed Loop  
% [t,yNL] = ode15s(@(t,y)roll_dynamics_4DOF_NL18(y,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,int,I,CL),tspan,x0); %Open Loop 

Draw

figure(4)
xlim([0 40])


%%
if SG==1
saveas(figure(1),'sim123a');saveas(figure(1),'sim123b.jpg');
saveas(figure(2),'sim123c');saveas(figure(2),'sim123d.jpg');
saveas(figure(3),'sim123e');saveas(figure(3),'sim123f.jpg');
saveas(figure(4),'sim123g');saveas(figure(4),'sim123h.jpg');
saveas(figure(6),'sim123i');saveas(figure(6),'sim123j.jpg');
saveas(figure(7),'sim123k');saveas(figure(7),'sim123l.jpg');
saveas(figure(8),'sim123m');saveas(figure(8),'sim123n.jpg');
saveas(figure(9),'sim123o');saveas(figure(9),'sim123p.jpg');
saveas(figure(5),'sim123q');saveas(figure(5),'sim123r.jpg');
save('sim_data')
end

%% Sensitity analysis routines for parameters zm, vcg
% % Evaluate linear controllability as function of zm
% for j=1:40
%     zm=j*.1;
% [AA,BB] = MIMO(x0,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,I,int,CL);
% B1=BB(:,1);
% CC1=W'*B1;
% CC1_j(j)=max(CC1);
% end
% plot(.1:.1:40*.1,CC1_j);
% ylabel('Controllability')
% xlabel('z_o')

% Evaluate linear controllability as function of vcg
% AA0=AA0(1:6,1:6);  %mass fixed 4DOF system
% [V0,E0,W0]=eig(AA0);
% for j=1:10
%     vcg=j*.1;
% [AA,BB] = MIMO(x0,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,I,int,CL);
% B1r=BB(1:6,2);
% CC1r=W0'*B1r;
% CC1r_j(j)=max(CC1r);
% end
% plot(CC1r_j);

% Controllability evaluation
% CC1_m=norm(CC);
% CC2=W'*B2;
% CC2_m=norm(CC);
% [V0,E0,W0]=eig(AA0); 
% CC1r=W0'*B1r;
% CC1r_m=norm(CC);
% [CC1_m CC2_m CC1r_m]
